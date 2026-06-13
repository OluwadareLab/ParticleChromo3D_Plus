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
