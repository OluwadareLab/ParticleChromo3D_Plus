pub mod helper;
pub mod swarm;

use log::{debug, info};
use rayon::prelude::*;
use swarm::{LossFunc, Swarm, loss_function, pearsonr, spearmanr};

#[derive(Clone, Copy)]
pub struct PsoParams {
    pub rand_range: f64,
    pub swarm_size: usize,
    pub threshold: f64,
    pub itt_count: usize,
    pub loss_func: LossFunc,
}

pub struct OptResult {
    pub pearson: f64,
    pub spearman: f64,
    pub cost: f64,
    pub itt_fin: usize,
    pub swarm_id: usize,
    pub swarm: Swarm,
    pub alpha_idx: usize,
}

pub fn alpha_grid(start: f64, end: f64, step: f64) -> Vec<f64> {
    let mut grid = vec![];
    let mut alpha = start;
    while alpha < end - 1e-9 {
        grid.push(alpha);
        alpha += step;
    }
    grid
}

fn one_move(swarm: &mut Swarm, target: &[f64], params: PsoParams) -> usize {
    let mut save_g_best_cost = f64::INFINITY;

    for i in 0..params.itt_count {
        if i % 1000 == 0
            && let Some(ref g) = swarm.g_best
        {
            let error = loss_function(target, &g.2, params.loss_func);
            debug!(
                "id: {} itt: {} Cost: {} Pearson: {} Spearman: {}",
                swarm.id,
                i,
                g.1,
                pearsonr(&g.2, target),
                spearmanr(&g.2, target),
            );

            if (save_g_best_cost - error).abs() >= params.threshold {
                save_g_best_cost = error;
            } else {
                return i;
            }
        }

        swarm.calc_vel(params.itt_count, i);
        swarm.update_pos(i);
        swarm.calc_cost();
    }

    params.itt_count - 1
}

pub fn optimize(
    contacts_with_dist: Vec<[f64; 4]>,
    point_count: usize,
    zero_ind: Vec<usize>,
    params: PsoParams,
    alpha_idx: usize,
) -> OptResult {
    let target: Vec<f64> = contacts_with_dist.iter().map(|c| c[3]).collect();

    let mut swarm = Swarm::new(
        contacts_with_dist,
        point_count,
        params.rand_range,
        params.swarm_size,
        zero_ind,
    );
    swarm.loss_func = params.loss_func;

    let itt_fin = one_move(&mut swarm, &target, params);

    let g = swarm.g_best.as_ref().unwrap();
    let pearson = pearsonr(&g.2, &target);
    let spearman = spearmanr(&g.2, &target);
    let cost = loss_function(&target, &g.2, params.loss_func);
    let swarm_id = swarm.id;

    OptResult {
        pearson,
        spearman,
        cost,
        itt_fin,
        swarm_id,
        swarm,
        alpha_idx,
    }
}

pub fn par_choice(
    file_ptr: &str,
    out_file_ptr: &str,
    alphas: &[f64],
    params: PsoParams,
) -> OptResult {
    let (contacts, point_map, zero_ind) = helper::read_data(file_ptr);
    let point_count = point_map.len();

    info!(
        "Running PSO over {} alpha values with {} threads",
        alphas.len(),
        rayon::current_num_threads()
    );

    let results: Vec<OptResult> = alphas
        .par_iter()
        .enumerate()
        .map(|(idx, &alpha)| {
            let contacts_4: Vec<[f64; 4]> = contacts
                .iter()
                .map(|c| [c[0], c[1], c[2], 1.0 / c[2].powf(alpha)])
                .collect();

            optimize(contacts_4, point_count, zero_ind.clone(), params, idx)
        })
        .collect();

    let best = results
        .into_iter()
        .max_by(|a, b| a.spearman.partial_cmp(&b.spearman).unwrap())
        .unwrap();

    let g = best.swarm.g_best.as_ref().unwrap();
    helper::write_output(out_file_ptr, &g.0);

    best
}
