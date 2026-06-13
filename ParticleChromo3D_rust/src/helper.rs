use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::Path;

/// Reads a space-delimited square contact matrix and returns:
/// - contact list as Vec<[f64; 3]> (i, j, freq)
/// - point map (original index -> compact index)
/// - zero_ind: indices in pdist output corresponding to zero-contact pairs
pub fn read_matrix_to_list(file_ptr: &str) -> (Vec<[f64; 3]>, Vec<usize>) {
    let content = fs::read_to_string(file_ptr).expect("Failed to read input file");

    // Parse rows
    let mut matrix: Vec<Vec<f64>> = content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            line.split_whitespace()
                .map(|v| v.parse::<f64>().unwrap_or(0.0))
                .collect()
        })
        .collect();

    // Delete all-zero rows
    matrix.retain(|row| row.iter().any(|&v| v != 0.0));

    if matrix.is_empty() {
        return (vec![], vec![]);
    }

    let n_rows = matrix.len();
    let n_cols = matrix[0].len();

    // Find all-zero columns
    let mut zero_cols: Vec<usize> = vec![];
    for col in 0..n_cols {
        if (0..n_rows).all(|row| matrix[row][col] == 0.0) {
            zero_cols.push(col);
        }
    }

    // Remove zero columns
    let keep_cols: Vec<usize> = (0..n_cols)
        .filter(|c| !zero_cols.contains(c))
        .collect();

    let matrix: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| keep_cols.iter().map(|&c| row[c]).collect())
        .collect();

    let n = matrix.len().min(matrix[0].len());

    let mut stop_dupe: HashSet<(usize, usize)> = HashSet::new();
    let mut contact_list: Vec<[f64; 3]> = vec![];
    let mut zero_ind: Vec<usize> = vec![];
    let mut count = 0usize;

    for i in 0..n {
        for j in 0..matrix[i].len() {
            if i != j && !stop_dupe.contains(&(i, j)) && !stop_dupe.contains(&(j, i)) {
                stop_dupe.insert((i, j));
                stop_dupe.insert((j, i));
                if matrix[i][j] > 0.0 {
                    contact_list.push([i as f64, j as f64, matrix[i][j]]);
                } else {
                    zero_ind.push(count);
                }
                count += 1;
            }
        }
    }

    let zero_ind = if zero_ind.is_empty() { vec![] } else { zero_ind };
    (contact_list, zero_ind)
}

/// Wraps read_matrix_to_list and remaps point indices to a compact range.
/// Returns (contact array, point_map, zero_ind).
pub fn read_data(file_ptr: &str) -> (Vec<[f64; 3]>, HashMap<usize, usize>, Vec<usize>) {
    let (mut contacts, zero_ind) = read_matrix_to_list(file_ptr);

    let mut point_set: HashSet<usize> = HashSet::new();
    for c in &contacts {
        point_set.insert(c[0] as usize);
        point_set.insert(c[1] as usize);
    }

    let mut point_map: HashMap<usize, usize> = HashMap::new();
    for (new_idx, orig) in point_set.into_iter().enumerate() {
        point_map.insert(orig, new_idx);
    }

    for c in contacts.iter_mut() {
        c[0] = *point_map.get(&(c[0] as usize)).unwrap() as f64;
        c[1] = *point_map.get(&(c[1] as usize)).unwrap() as f64;
    }

    (contacts, point_map, zero_ind)
}

/// Strips excess whitespace from a file, writes to <file>.stripped
pub fn strip_file(in_file: &str) -> String {
    let content = fs::read_to_string(in_file).expect("Failed to read file for stripping");
    let cleaned: Vec<String> = content
        .lines()
        .map(|line| line.split_whitespace().collect::<Vec<_>>().join(" "))
        .collect();
    let fout = format!("{}.stripped", in_file);
    let mut f = fs::File::create(&fout).expect("Failed to create stripped file");
    writeln!(f, "{}", cleaned.join("\n")).unwrap();
    fout
}

/// Scales xyz array to [min_val, max_val]
pub fn scale_arr(xyz: &mut Vec<[f64; 3]>, min_val: f64, max_val: f64) {
    let flat_min = xyz.iter().flat_map(|p| p.iter().copied()).fold(f64::INFINITY, f64::min);
    let flat_max = xyz.iter().flat_map(|p| p.iter().copied()).fold(f64::NEG_INFINITY, f64::max);
    let old_range = flat_max - flat_min;
    let new_range = max_val - min_val;
    for p in xyz.iter_mut() {
        for v in p.iter_mut() {
            *v = ((*v - flat_min) * new_range / old_range) + min_val;
        }
    }
}

/// Writes xyz coordinates as a PDB file
pub fn write_pdb(positions: &[[f64; 3]], pdb_file: &str) {
    let out_dir = Path::new("./out");
    fs::create_dir_all(out_dir).ok();

    let mut f = fs::File::create(pdb_file).expect("Failed to create PDB file");
    writeln!(f).unwrap();

    let bin_num = positions.len();
    for (idx, pos) in positions.iter().enumerate() {
        let i = idx + 1;
        let col1 = "ATOM";
        let col2 = format!("{:>5}", i);
        let col3 = "CA MET";
        let col4 = format!("{:<6}", format!("B{}", i));
        let col5 = format!("{:>8.3}", pos[0]);
        let col6 = format!("{:>8.3}", pos[1]);
        let col7 = format!("{:>8.3}", pos[2]);
        let col8 = "0.20 10.00";
        writeln!(f, "{}  {}   {} {}   {}{}{}  {}", col1, col2, col3, col4, col5, col6, col7, col8).unwrap();
    }

    for i in 1..=bin_num {
        let j = i + 1;
        let line = format!("CONECT{:>5}{:>5}", i, j);
        writeln!(f, "{}", line).unwrap();
    }
    writeln!(f, "END").unwrap();
}

/// Writes xyz to a PDB (scales first)
pub fn write_output(file_ptr: &str, xyz: &[[f64; 3]]) {
    let mut xyz_owned = xyz.to_vec();
    scale_arr(&mut xyz_owned, -10.0, 10.0);
    write_pdb(&xyz_owned, &format!("{}.pdb", file_ptr));
}

/// Writes a summary log file
pub fn write_log(outfile: &str, in_file: &str, best_alpha: f64, rmse: f64, best_spearman: f64, best_pearson: f64) {
    let content = format!(
        "Input file: {}\nConvert factor:: {}\nBest cost  : {}\nBest Spearman correlation Dist vs. Reconstructed Dist  : {}\nBest Pearson correlation Dist vs. Reconstructed Dist  : {}\n",
        in_file, best_alpha, rmse, best_spearman, best_pearson
    );
    fs::write(outfile, content).expect("Failed to write log");
}
