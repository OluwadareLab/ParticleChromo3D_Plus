use clap::Parser;
use log::info;
use particle_chromo3d::swarm::LossFunc;
use particle_chromo3d::{PsoParams, alpha_grid, helper, par_choice};
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

    let alphas = alpha_grid(0.1, 2.0, 0.1);

    let params = PsoParams {
        rand_range: args.rand_range,
        swarm_size: args.swarm_size,
        threshold: args.threshold,
        itt_count: args.itt_count,
        loss_func,
    };

    let best = par_choice(&stripped, &out_file_ptr, &alphas, params);

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
