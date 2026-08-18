use particle_chromo3d::swarm::{LossFunc, pdist, spearmanr};
use particle_chromo3d::{PsoParams, alpha_grid, par_choice};
use std::fs;
use tempfile::tempdir;

fn helix_points(n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            let t = i as f64 * 0.6;
            [t.cos() * 5.0, t.sin() * 5.0, i as f64 * 1.2]
        })
        .collect()
}

fn contact_matrix(points: &[[f64; 3]]) -> String {
    let n = points.len();
    let mut rows = Vec::with_capacity(n);
    for i in 0..n {
        let cells: Vec<String> = (0..n)
            .map(|j| {
                if i == j {
                    "0".to_string()
                } else {
                    let dx = points[i][0] - points[j][0];
                    let dy = points[i][1] - points[j][1];
                    let dz = points[i][2] - points[j][2];
                    format!("{:.6}", 1.0 / (dx * dx + dy * dy + dz * dz).sqrt())
                }
            })
            .collect();
        rows.push(cells.join(" "));
    }
    rows.join("\n") + "\n"
}

fn atom_coordinates(pdb: &str) -> Vec<[f64; 3]> {
    pdb.lines()
        .filter(|l| l.starts_with("ATOM"))
        .map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            let n = t.len();
            [
                t[n - 5].parse().unwrap(),
                t[n - 4].parse().unwrap(),
                t[n - 3].parse().unwrap(),
            ]
        })
        .collect()
}

#[test]
fn par_choice_reconstructs_a_known_helix_and_writes_a_valid_pdb() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("helix.txt");
    let out_base = dir.path().join("chr");
    let truth = helix_points(12);
    fs::write(&input, contact_matrix(&truth)).unwrap();

    let params = PsoParams {
        rand_range: 1.0,
        swarm_size: 5,
        threshold: 1e-6,
        itt_count: 2000,
        loss_func: LossFunc::Rmse,
    };
    let alphas = alpha_grid(0.1, 1.0, 0.1);

    let best = par_choice(
        input.to_str().unwrap(),
        out_base.to_str().unwrap(),
        &alphas,
        params,
    );

    assert!(
        best.spearman > 0.9,
        "reconstruction should track the true structure, got spearman {}",
        best.spearman
    );
    assert!(best.alpha_idx < alphas.len());
    assert!(best.cost.is_finite());

    let pdb = fs::read_to_string(format!("{}.pdb", out_base.to_str().unwrap())).unwrap();
    let coords = atom_coordinates(&pdb);

    assert_eq!(coords.len(), truth.len());
    assert_eq!(pdb.lines().filter(|l| l.starts_with("CONECT")).count(), 12);
    assert_eq!(pdb.trim_end().lines().last().unwrap(), "END");

    for bead in &coords {
        for &v in bead {
            assert!((-10.0..=10.0).contains(&v), "coordinate {v} left the box");
        }
    }

    let recovered = spearmanr(&pdist(&coords), &pdist(&truth));
    assert!(
        recovered > 0.9,
        "pdb geometry should match the true helix, got spearman {recovered}"
    );
}
