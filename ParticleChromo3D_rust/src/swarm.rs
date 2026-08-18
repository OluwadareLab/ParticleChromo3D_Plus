use rand::Rng;
use std::sync::atomic::{AtomicUsize, Ordering};

static SWARM_ID_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// Loss function choices matching the Python implementation
#[derive(Clone, Copy, Debug)]
pub enum LossFunc {
    Sse = 0,
    Mse = 1,
    Rmse = 2,
    Huber = 3,
}

impl LossFunc {
    pub fn from_int(n: i32) -> LossFunc {
        match n {
            0 => LossFunc::Sse,
            1 => LossFunc::Mse,
            2 => LossFunc::Rmse,
            3 => LossFunc::Huber,
            _ => LossFunc::Rmse,
        }
    }
}

/// Compute pairwise euclidean distances (upper triangle, row-major order)
pub fn pdist(positions: &[[f64; 3]]) -> Vec<f64> {
    let n = positions.len();
    let mut dists = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            dists.push((dx * dx + dy * dy + dz * dz).sqrt());
        }
    }
    dists
}

/// Pearson correlation coefficient
pub fn pearsonr(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let num: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let dx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum::<f64>().sqrt();
    let dy: f64 = y.iter().map(|b| (b - my).powi(2)).sum::<f64>().sqrt();
    if dx == 0.0 || dy == 0.0 { return 0.0; }
    num / (dx * dy)
}

/// Spearman rank correlation coefficient
pub fn spearmanr(x: &[f64], y: &[f64]) -> f64 {
    let rank = |v: &[f64]| -> Vec<f64> {
        let mut indexed: Vec<(usize, f64)> = v.iter().copied().enumerate().collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let mut ranks = vec![0.0f64; v.len()];
        let mut i = 0;
        while i < indexed.len() {
            let mut j = i;
            while j + 1 < indexed.len() && indexed[j + 1].1 == indexed[i].1 {
                j += 1;
            }
            let avg_rank = (i + j) as f64 / 2.0 + 1.0;
            for k in i..=j {
                ranks[indexed[k].0] = avg_rank;
            }
            i = j + 1;
        }
        ranks
    };
    pearsonr(&rank(x), &rank(y))
}

/// Compute loss between predicted distances and target distances
pub fn loss_function(target: &[f64], predicted: &[f64], func: LossFunc) -> f64 {
    let n = target.len() as f64;
    match func {
        LossFunc::Rmse => {
            let sum: f64 = target.iter().zip(predicted).map(|(t, p)| (p - t).powi(2)).sum();
            (sum / n).sqrt()
        }
        LossFunc::Sse => {
            target.iter().zip(predicted).map(|(t, p)| (p - t).powi(2)).sum()
        }
        LossFunc::Mse => {
            let sum: f64 = target.iter().zip(predicted).map(|(t, p)| (p - t).powi(2)).sum();
            sum / n
        }
        LossFunc::Huber => {
            let alpha = 0.5f64;
            target.iter().zip(predicted).map(|(t, p)| {
                let diff = (p - t).abs();
                if diff < alpha { 0.5 * (p - t).powi(2) }
                else { alpha * (diff - 0.5 * alpha) }
            }).sum()
        }
    }
}

/// A particle is a flat vec of xyz for each bead: shape [point_count][3]
type Particle = Vec<[f64; 3]>;

pub struct Swarm {
    pub id: usize,
    pub pc: usize,
    pub rand_max: f64,
    pub rand_min: f64,
    pub loss_func: LossFunc,

    /// Global best: (positions, cost, distances)
    pub g_best: Option<(Vec<[f64; 3]>, f64, Vec<f64>)>,

    /// ref_contacts: Vec<[i, j, if, target_dist]> — 4-column
    pub ref_contacts: Vec<[f64; 4]>,
    pub zero_ind: Vec<usize>,

    /// pos[particle][bead][xyz]
    pub pos: Vec<Particle>,
    pub pos_best: Vec<Particle>,
    pub cost_best: Vec<f64>,
    pub vel: Vec<Particle>,
    pub cost: Vec<f64>,
    pub dist: Vec<Vec<f64>>,
    pub loc_op_count: Vec<f64>,
}

impl Swarm {
    pub fn new(
        ref_contacts: Vec<[f64; 4]>,
        point_count: usize,
        rand_val: f64,
        swarm_size: usize,
        zero_ind: Vec<usize>,
    ) -> Self {
        let id = SWARM_ID_COUNTER.fetch_add(1, Ordering::SeqCst) + 1;
        let rand_max = rand_val;
        let rand_min = -rand_val;

        let mut rng = rand::thread_rng();
        let pos: Vec<Particle> = (0..swarm_size)
            .map(|_| Self::rand_cur_static(&mut rng, point_count, rand_min, rand_max))
            .collect();

        let n_dists = ref_contacts.len();
        let pos_best = pos.clone();
        let cost_best = vec![f64::INFINITY; swarm_size];
        let vel: Vec<Particle> = (0..swarm_size)
            .map(|_| vec![[0.0; 3]; point_count])
            .collect();
        let cost = vec![f64::INFINITY; swarm_size];
        let dist = vec![vec![0.0; n_dists]; swarm_size];
        let loc_op_count = vec![0.0; swarm_size];

        let mut swarm = Swarm {
            id,
            pc: point_count,
            rand_max,
            rand_min,
            loss_func: LossFunc::Rmse,
            g_best: None,
            ref_contacts,
            zero_ind,
            pos,
            pos_best,
            cost_best,
            vel,
            cost,
            dist,
            loc_op_count,
        };

        swarm.calc_cost();
        swarm
    }

    fn rand_cur_static(rng: &mut impl Rng, pc: usize, min: f64, max: f64) -> Particle {
        (0..pc).map(|_| {
            [
                rng.gen_range(min..=max),
                rng.gen_range(min..=max),
                rng.gen_range(min..=max),
            ]
        }).collect()
    }

    fn rand_cur(&self, rng: &mut impl Rng) -> Particle {
        Self::rand_cur_static(rng, self.pc, self.rand_min, self.rand_max)
    }

    fn rand_shift(rng: &mut impl Rng, copy_pos: &Particle, cut_size: usize, threshold: f64) -> (Particle, Vec<bool>) {
        let mut temp = copy_pos.clone();
        let n = temp.len();
        let mut mask = vec![false; n];

        // create a boolean mask: cut_size falses, rest trues, then shuffle
        for i in cut_size..n {
            mask[i] = true;
        }
        // Fisher-Yates shuffle on mask
        for i in (1..n).rev() {
            let j = rng.gen_range(0..=i);
            mask.swap(i, j);
        }

        for i in 0..n {
            if mask[i] {
                for k in 0..3 {
                    temp[i][k] += rng.gen_range(-threshold..=threshold);
                }
            }
        }
        (temp, mask)
    }

    pub fn calc_dist(&mut self) {
        for p in 0..self.pos.len() {
            let full = pdist(&self.pos[p]);
            if self.zero_ind.is_empty() {
                self.dist[p] = full;
            } else {
                self.dist[p] = full.into_iter().enumerate()
                    .filter(|(i, _)| !self.zero_ind.contains(i))
                    .map(|(_, v)| v)
                    .collect();
            }
        }
    }

    fn compute_cost_for_particle(&self, p: usize) -> f64 {
        let target: Vec<f64> = self.ref_contacts.iter().map(|c| c[3]).collect();
        let predicted = &self.dist[p];
        match self.loss_func {
            LossFunc::Rmse => {
                let sum: f64 = target.iter().zip(predicted).map(|(t, d)| (d - t).powi(2)).sum();
                sum.sqrt()
            }
            LossFunc::Sse => {
                target.iter().zip(predicted).map(|(t, d)| (d - t).powi(2)).sum()
            }
            LossFunc::Mse => {
                let n = target.len() as f64;
                let sum: f64 = target.iter().zip(predicted).map(|(t, d)| (d - t).powi(2)).sum();
                (1.0 / self.pc as f64) * sum / n
            }
            LossFunc::Huber => {
                let delta = 0.1f64;
                target.iter().zip(predicted).map(|(t, d)| {
                    let diff = (d - t).abs();
                    if diff < delta { 0.5 * (d - t).powi(2) }
                    else { delta * (diff - 0.5 * delta) }
                }).sum()
            }
        }
    }

    pub fn calc_cost(&mut self) {
        self.calc_dist();
        let n = self.pos.len();
        let new_costs: Vec<f64> = (0..n).map(|p| self.compute_cost_for_particle(p)).collect();
        self.update_cost(new_costs);
    }

    fn update_cost(&mut self, new_cost: Vec<f64>) {
        for p in 0..self.pos.len() {
            if new_cost[p] > self.cost[p] {
                self.loc_op_count[p] += 1.0;
            }
            self.cost[p] = new_cost[p];

            if self.cost[p] < self.cost_best[p] {
                self.pos_best[p] = self.pos[p].clone();
                self.cost_best[p] = self.cost[p];
            }
        }

        let best_p = self.cost.iter().enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        let best_cost = self.cost[best_p];
        let is_new_best = self.g_best.as_ref().map_or(true, |g| best_cost < g.1);
        if is_new_best {
            self.g_best = Some((
                self.pos[best_p].clone(),
                best_cost,
                self.dist[best_p].clone(),
            ));
        }
    }

    fn calc_const(itt_max: f64, k: f64, w_min: f64, w_max: f64) -> f64 {
        let k = k.min(itt_max);
        ((w_max - w_min) * ((itt_max - k) / itt_max)) + w_min
    }

    pub fn calc_vel(&mut self, _itt_max: usize, _itt: usize) {
        let weight = 0.5f64;
        let con_p = 0.3f64;
        let con_g = 2.5f64;
        let g_best_pos = self.g_best.as_ref().unwrap().0.clone();
        let mut rng = rand::thread_rng();

        for p in 0..self.pos.len() {
            for b in 0..self.pc {
                for k in 0..3 {
                    let ran_p: f64 = rng.r#gen();
                    let ran_g: f64 = rng.r#gen();
                    self.vel[p][b][k] = weight * self.vel[p][b][k]
                        + con_p * ran_p * (self.pos_best[p][b][k] - self.pos[p][b][k])
                        + con_g * ran_g * (g_best_pos[b][k] - self.pos[p][b][k]);
                }
            }
        }
    }

    pub fn update_pos(&mut self, itt: usize) {
        let mut rng = rand::thread_rng();
        let cut_size = rng.gen_range(1..self.pc.saturating_sub(1).max(2));
        let thresh = if itt > 500 { (1.0 / itt as f64) * 100.0 } else { 1.0 };
        let const_val = Self::calc_const(10000.0, itt as f64, 5.0, 15.0);

        let changed: Vec<usize> = (0..self.pos.len())
            .filter(|&p| self.loc_op_count[p] > const_val)
            .collect();

        for &p in &changed {
            if itt < 1000 {
                self.pos[p] = self.rand_cur(&mut rng);
            } else {
                let (new_pos, mask) = Self::rand_shift(&mut rng, &self.pos[p].clone(), cut_size, thresh);
                self.pos[p] = new_pos;
                for b in 0..self.pc {
                    if mask[b] {
                        self.vel[p][b] = [0.0; 3];
                    }
                }
            }
            self.vel[p] = vec![[0.0; 3]; self.pc];
            self.loc_op_count[p] = -1.0;
        }

        let not_changed: Vec<usize> = (0..self.pos.len())
            .filter(|p| !changed.contains(p))
            .collect();

        for &p in &not_changed {
            for b in 0..self.pc {
                for k in 0..3 {
                    self.pos[p][b][k] += self.vel[p][b][k];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-9,
            "expected {expected}, got {actual}"
        );
    }

    fn three_point_swarm(swarm_size: usize) -> Swarm {
        Swarm::new(
            vec![
                [0.0, 1.0, 1.0, 1.0],
                [0.0, 2.0, 1.0, 2.0],
                [1.0, 2.0, 1.0, 1.5],
            ],
            3,
            1.0,
            swarm_size,
            vec![],
        )
    }

    fn two_point_swarm() -> Swarm {
        Swarm::new(vec![[0.0, 1.0, 1.0, 5.0]], 2, 1.0, 1, vec![])
    }

    #[test]
    fn loss_func_from_int_maps_known_codes_and_falls_back_to_rmse() {
        assert_eq!(LossFunc::from_int(0) as i32, LossFunc::Sse as i32);
        assert_eq!(LossFunc::from_int(1) as i32, LossFunc::Mse as i32);
        assert_eq!(LossFunc::from_int(2) as i32, LossFunc::Rmse as i32);
        assert_eq!(LossFunc::from_int(3) as i32, LossFunc::Huber as i32);
        assert_eq!(LossFunc::from_int(-7) as i32, LossFunc::Rmse as i32);
        assert_eq!(LossFunc::from_int(99) as i32, LossFunc::Rmse as i32);
    }

    #[test]
    fn pdist_returns_upper_triangle_distances_in_row_major_order() {
        let dists = pdist(&[[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 1.0]]);

        assert_eq!(dists.len(), 3);
        assert_close(dists[0], 5.0);
        assert_close(dists[1], 1.0);
        assert_close(dists[2], 26f64.sqrt());
    }

    #[test]
    fn pdist_of_a_single_point_is_empty() {
        assert!(pdist(&[[1.0, 2.0, 3.0]]).is_empty());
    }

    #[test]
    fn pearsonr_detects_perfect_positive_and_negative_linear_relationships() {
        assert_close(pearsonr(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]), 1.0);
        assert_close(pearsonr(&[1.0, 2.0, 3.0], &[3.0, 2.0, 1.0]), -1.0);
    }

    #[test]
    fn pearsonr_returns_zero_when_a_series_has_no_variance() {
        assert_close(pearsonr(&[1.0, 1.0, 1.0], &[1.0, 2.0, 3.0]), 0.0);
        assert_close(pearsonr(&[1.0, 2.0, 3.0], &[4.0, 4.0, 4.0]), 0.0);
    }

    #[test]
    fn spearmanr_is_one_for_any_monotonic_relationship() {
        assert_close(spearmanr(&[1.0, 2.0, 3.0], &[1.0, 10.0, 1000.0]), 1.0);
        assert_close(spearmanr(&[1.0, 2.0, 3.0], &[1000.0, 10.0, 1.0]), -1.0);
    }

    #[test]
    fn spearmanr_assigns_tied_values_their_average_rank() {
        assert_close(
            spearmanr(&[10.0, 20.0, 20.0, 30.0], &[1.0, 2.0, 3.0, 4.0]),
            0.9f64.sqrt(),
        );
    }

    #[test]
    fn loss_function_is_zero_for_an_exact_match() {
        let v = [1.0, 2.0, 3.0];
        for func in [LossFunc::Sse, LossFunc::Mse, LossFunc::Rmse, LossFunc::Huber] {
            assert_close(loss_function(&v, &v, func), 0.0);
        }
    }

    #[test]
    fn loss_function_computes_each_squared_error_variant() {
        let target = [0.0, 0.0];
        let predicted = [1.0, 3.0];

        assert_close(loss_function(&target, &predicted, LossFunc::Sse), 10.0);
        assert_close(loss_function(&target, &predicted, LossFunc::Mse), 5.0);
        assert_close(loss_function(&target, &predicted, LossFunc::Rmse), 5f64.sqrt());
    }

    #[test]
    fn loss_function_huber_switches_from_quadratic_to_linear_at_alpha() {
        assert_close(loss_function(&[0.0], &[0.2], LossFunc::Huber), 0.02);
        assert_close(loss_function(&[0.0], &[3.0], LossFunc::Huber), 1.375);
        assert_close(loss_function(&[0.0, 0.0], &[1.0, 3.0], LossFunc::Huber), 1.75);
    }

    #[test]
    fn new_swarm_allocates_consistent_shapes_for_every_particle() {
        let swarm = three_point_swarm(4);

        assert_eq!(swarm.pos.len(), 4);
        assert_eq!(swarm.pos_best.len(), 4);
        assert_eq!(swarm.vel.len(), 4);
        for p in 0..4 {
            assert_eq!(swarm.pos[p].len(), 3);
            assert_eq!(swarm.vel[p], vec![[0.0; 3]; 3]);
            assert_eq!(swarm.dist[p].len(), 3);
        }
    }

    #[test]
    fn new_swarm_seeds_positions_inside_the_requested_range() {
        let swarm = three_point_swarm(8);

        assert_close(swarm.rand_min, -1.0);
        assert_close(swarm.rand_max, 1.0);
        for particle in &swarm.pos {
            for bead in particle {
                for &v in bead {
                    assert!((-1.0..=1.0).contains(&v), "coordinate {v} out of range");
                }
            }
        }
    }

    #[test]
    fn new_swarm_scores_every_particle_and_records_the_global_best() {
        let swarm = three_point_swarm(5);

        assert!(swarm.cost.iter().all(|c| c.is_finite()));
        assert_eq!(swarm.cost_best, swarm.cost);
        assert!(swarm.loc_op_count.iter().all(|&c| c == 0.0));

        let (_, best_cost, best_dist) = swarm.g_best.as_ref().unwrap();
        let min_cost = swarm.cost.iter().copied().fold(f64::INFINITY, f64::min);
        assert_close(*best_cost, min_cost);
        assert_eq!(best_dist.len(), 3);
    }

    #[test]
    fn swarm_ids_are_unique() {
        let a = three_point_swarm(2);
        let b = three_point_swarm(2);

        assert_ne!(a.id, b.id);
    }

    #[test]
    fn calc_dist_keeps_every_pair_when_no_zero_indices_are_given() {
        let mut swarm = three_point_swarm(1);
        swarm.pos[0] = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 1.0]];

        swarm.calc_dist();

        assert_eq!(swarm.dist[0].len(), 3);
        assert_close(swarm.dist[0][0], 5.0);
        assert_close(swarm.dist[0][1], 1.0);
        assert_close(swarm.dist[0][2], 26f64.sqrt());
    }

    #[test]
    fn calc_dist_drops_the_pairs_listed_in_zero_ind() {
        let mut swarm = Swarm::new(
            vec![[0.0, 1.0, 1.0, 1.0], [1.0, 2.0, 1.0, 1.0]],
            3,
            1.0,
            1,
            vec![1],
        );
        swarm.pos[0] = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 1.0]];

        swarm.calc_dist();

        assert_eq!(swarm.dist[0].len(), 2);
        assert_close(swarm.dist[0][0], 5.0);
        assert_close(swarm.dist[0][1], 26f64.sqrt());
    }

    #[test]
    fn calc_cost_reaches_zero_when_a_particle_matches_the_target_distances() {
        let mut swarm = two_point_swarm();
        swarm.pos[0] = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];

        swarm.calc_cost();

        assert_close(swarm.cost[0], 0.0);
        assert_close(swarm.g_best.as_ref().unwrap().1, 0.0);
    }

    #[test]
    fn calc_cost_counts_regressions_without_losing_the_personal_best() {
        let mut swarm = two_point_swarm();
        let solution = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        swarm.pos[0] = solution.clone();
        swarm.calc_cost();

        swarm.pos[0] = vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        swarm.calc_cost();

        assert_close(swarm.cost[0], 95.0);
        assert_close(swarm.loc_op_count[0], 1.0);
        assert_close(swarm.cost_best[0], 0.0);
        assert_eq!(swarm.pos_best[0], solution);
        assert_close(swarm.g_best.as_ref().unwrap().1, 0.0);
    }

    #[test]
    fn calc_const_decays_linearly_from_w_max_to_w_min() {
        assert_close(Swarm::calc_const(10000.0, 0.0, 5.0, 15.0), 15.0);
        assert_close(Swarm::calc_const(10000.0, 5000.0, 5.0, 15.0), 10.0);
        assert_close(Swarm::calc_const(10000.0, 10000.0, 5.0, 15.0), 5.0);
    }

    #[test]
    fn calc_const_clamps_iterations_beyond_the_maximum() {
        assert_close(Swarm::calc_const(10000.0, 99999.0, 5.0, 15.0), 5.0);
    }

    #[test]
    fn calc_vel_leaves_velocity_at_zero_when_every_particle_sits_on_the_best() {
        let mut swarm = three_point_swarm(3);
        let shared = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]];
        for p in 0..swarm.pos.len() {
            swarm.pos[p] = shared.clone();
            swarm.pos_best[p] = shared.clone();
        }
        swarm.g_best = Some((shared.clone(), 0.0, pdist(&shared)));

        swarm.calc_vel(100, 0);

        for particle in &swarm.vel {
            assert_eq!(particle, &vec![[0.0; 3]; 3]);
        }
    }

    #[test]
    fn calc_vel_pulls_a_lagging_particle_toward_the_global_best() {
        let mut swarm = three_point_swarm(1);
        let behind = vec![[0.0, 0.0, 0.0]; 3];
        let ahead = vec![[10.0, 10.0, 10.0]; 3];
        swarm.pos[0] = behind.clone();
        swarm.pos_best[0] = behind;
        swarm.g_best = Some((ahead.clone(), 0.0, pdist(&ahead)));

        swarm.calc_vel(100, 0);

        for bead in &swarm.vel[0] {
            for &v in bead {
                assert!(v > 0.0, "expected a positive pull toward the best, got {v}");
                assert!(v <= 25.0 + EPS, "velocity {v} exceeds con_g * distance");
            }
        }
    }

    #[test]
    fn update_pos_adds_velocity_to_particles_that_are_still_improving() {
        let mut swarm = three_point_swarm(2);
        let start = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
        for p in 0..swarm.pos.len() {
            swarm.pos[p] = start.clone();
            swarm.vel[p] = vec![[0.5, -0.5, 1.0]; 3];
        }

        swarm.update_pos(0);

        for p in 0..swarm.pos.len() {
            for b in 0..3 {
                assert_close(swarm.pos[p][b][0], start[b][0] + 0.5);
                assert_close(swarm.pos[p][b][1], start[b][1] - 0.5);
                assert_close(swarm.pos[p][b][2], start[b][2] + 1.0);
            }
        }
    }

    #[test]
    fn update_pos_restarts_particles_that_have_stagnated() {
        let mut swarm = three_point_swarm(2);
        swarm.vel[0] = vec![[9.0, 9.0, 9.0]; 3];
        swarm.loc_op_count[0] = 1000.0;

        swarm.update_pos(0);

        assert_close(swarm.loc_op_count[0], -1.0);
        assert_eq!(swarm.vel[0], vec![[0.0; 3]; 3]);
        for bead in &swarm.pos[0] {
            for &v in bead {
                assert!((-1.0..=1.0).contains(&v), "restarted coordinate {v} out of range");
            }
        }
        assert_close(swarm.loc_op_count[1], 0.0);
    }

    #[test]
    fn rand_shift_perturbs_exactly_the_masked_beads_within_the_threshold() {
        let mut rng = rand::thread_rng();
        let original: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, -(i as f64)]).collect();
        let cut_size = 4;
        let threshold = 0.25;

        let (shifted, mask) = Swarm::rand_shift(&mut rng, &original, cut_size, threshold);

        assert_eq!(shifted.len(), original.len());
        assert_eq!(mask.iter().filter(|&&m| m).count(), original.len() - cut_size);
        for i in 0..original.len() {
            for k in 0..3 {
                let delta = shifted[i][k] - original[i][k];
                if mask[i] {
                    assert!(delta.abs() <= threshold + EPS, "delta {delta} exceeds threshold");
                } else {
                    assert_close(delta, 0.0);
                }
            }
        }
    }
}
