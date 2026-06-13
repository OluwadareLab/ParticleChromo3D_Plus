mod helper;
mod swarm;

use clap::Parser;
use log::{debug, info};
use rayon::prelude::*;
use swarm::{loss_function, pearsonr, spearmanr, LossFunc, Swarm};
use uuid::Uuid;

#[derive(Parser, Debug)]
#[command(name = "ParticleChromo3D")]
struct Args {
    /// Matrix of contacts (input file)
    infile: String,

    /// Number of particles in the system
    #[arg(short = 's', long = "swarmSize", default_value_t = 5)]
    swarm_size: usize,

    /// Maximum iterations before stop
    #[arg(short = 'i', long = "ittCount", default_value_t = 30000)]
    itt_count: usize,

    /// Error threshold before stopping
    #[arg(short = 't', long = "threshold", default_value_t = 0.000001)]
    threshold: f64,

    /// Range of x,y,z starting coords (uniform [-randRange, randRange])
    #[arg(short = 'r', long = "randRange", default_value_t = 1.0)]
    rand_range: f64,

    /// Output PDB filename
    #[arg(short = 'o', long = "outfile", default_value = "./out/chr.pdb")]
    outfile: String,

    /// Loss function (0=SSE, 1=MSE, 2=RMSE, 3=Huber)
    #[arg(short = 'l', long = "lossFunction", default_value_t = 2)]
    loss_function: i32,

    /// Log level (error, warn, info, debug, trace)
    #[arg(long = "logLevel", default_value = "info")]
    log_level: String,
}

struct OptResult {
    pearson: f64,
    spearman: f64,
    cost: f64,
    #[allow(dead_code)]
    itt_fin: usize,
    #[allow(dead_code)]
    swarm_id: usize,
    swarm: Swarm,
    alpha_idx: usize,
}

fn one_move(
    itt_count: usize,
    swarm: &mut Swarm,
    target: &[f64],
    threshold: f64,
    loss_func: LossFunc,
) -> usize {
    let mut save_g_best_cost = f64::INFINITY;

    for i in 0..itt_count {
        if i % 1000 == 0 {
            if let Some(ref g) = swarm.g_best {
                let error = loss_function(target, &g.2, loss_func);
                debug!(
                    "id: {} itt: {} Cost: {} Pearson: {} Spearman: {}",
                    swarm.id,
                    i,
                    g.1,
                    pearsonr(&g.2, target),
                    spearmanr(&g.2, target),
                );

                if (save_g_best_cost - error).abs() >= threshold {
                    save_g_best_cost = error;
                } else {
                    return i;
                }
            }
        }

        swarm.calc_vel(itt_count, i);
        swarm.update_pos(i);
        swarm.calc_cost();
    }

    itt_count - 1
}

fn optimize(
    contacts_with_dist: Vec<[f64; 4]>,
    point_count: usize,
    zero_ind: Vec<usize>,
    rand_range: f64,
    swarm_size: usize,
    threshold: f64,
    itt_count: usize,
    loss_func: LossFunc,
    alpha_idx: usize,
) -> OptResult {
    let target: Vec<f64> = contacts_with_dist.iter().map(|c| c[3]).collect();

    let mut swarm = Swarm::new(contacts_with_dist, point_count, rand_range, swarm_size, zero_ind);
    swarm.loss_func = loss_func;

    let itt_fin = one_move(itt_count, &mut swarm, &target, threshold, loss_func);

    let g = swarm.g_best.as_ref().unwrap();
    let pearson = pearsonr(&g.2, &target);
    let spearman = spearmanr(&g.2, &target);
    let cost = loss_function(&target, &g.2, loss_func);
    let swarm_id = swarm.id;

    OptResult { pearson, spearman, cost, itt_fin, swarm_id, swarm, alpha_idx }
}

fn par_choice(
    file_ptr: &str,
    out_file_ptr: &str,
    alpha_start: f64,
    alpha_end: f64,
    alpha_step: f64,
    rand_range: f64,
    swarm_size: usize,
    threshold: f64,
    itt_count: usize,
    loss_func: LossFunc,
) -> OptResult {
    let (contacts, point_map, zero_ind) = helper::read_data(file_ptr);
    let point_count = point_map.len();

    let alphas: Vec<f64> = {
        let mut v = vec![];
        let mut a = alpha_start;
        while a < alpha_end - 1e-9 {
            v.push(a);
            a += alpha_step;
        }
        v
    };

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

            optimize(
                contacts_4,
                point_count,
                zero_ind.clone(),
                rand_range,
                swarm_size,
                threshold,
                itt_count,
                loss_func,
                idx,
            )
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

fn main() {
    let args = Args::parse();

    let log_level = match args.log_level.to_uppercase().as_str() {
        "DEBUG" => log::LevelFilter::Debug,
        "WARNING" | "WARN" => log::LevelFilter::Warn,
        "ERROR" => log::LevelFilter::Error,
        "CRITICAL" => log::LevelFilter::Error,
        _ => log::LevelFilter::Info,
    };
    env_logger::Builder::new().filter_level(log_level).init();

    let out_file_ptr = format!("{}{}", args.outfile, Uuid::new_v4());
    let loss_func = LossFunc::from_int(args.loss_function);

    info!(
        "Starting ParticleChromo3D with log level: {}",
        args.log_level
    );
    info!("Processing file: {}", args.infile);

    let stripped = helper::strip_file(&args.infile);

    let alpha_start = 0.1f64;
    let alpha_end = 2.0f64;
    let alpha_step = 0.1f64;

    let best = par_choice(
        &stripped,
        &out_file_ptr,
        alpha_start,
        alpha_end,
        alpha_step,
        args.rand_range,
        args.swarm_size,
        args.threshold,
        args.itt_count,
        loss_func,
    );

    let alphas: Vec<f64> = {
        let mut v = vec![];
        let mut a = alpha_start;
        while a < alpha_end - 1e-9 {
            v.push(a);
            a += alpha_step;
        }
        v
    };
    let best_alpha = alphas[best.alpha_idx];

    info!("Input file: {}", args.infile);
    info!("Convert factor: {}", best_alpha);
    info!("Best cost: {}", best.cost);
    info!(
        "Best Spearman correlation Dist vs. Reconstructed Dist: {}",
        best.spearman
    );
    info!(
        "Best Pearson correlation Dist vs. Reconstructed Dist: {}",
        best.pearson
    );

    helper::write_log(
        &format!("{}.log", out_file_ptr),
        &args.infile,
        best_alpha,
        best.cost,
        best.spearman,
        best.pearson,
    );
}
