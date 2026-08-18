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

    let n_cols = matrix[0].len();

    // Find all-zero columns
    let zero_cols: Vec<usize> = (0..n_cols)
        .filter(|&col| matrix.iter().all(|row| row[col] == 0.0))
        .collect();

    // Remove zero columns
    let keep_cols: Vec<usize> = (0..n_cols).filter(|c| !zero_cols.contains(c)).collect();

    let matrix: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| keep_cols.iter().map(|&c| row[c]).collect())
        .collect();

    let n = matrix.len().min(matrix[0].len());

    let mut stop_dupe: HashSet<(usize, usize)> = HashSet::new();
    let mut contact_list: Vec<[f64; 3]> = vec![];
    let mut zero_ind: Vec<usize> = vec![];
    let mut count = 0usize;

    for (i, row) in matrix.iter().enumerate().take(n) {
        for (j, &value) in row.iter().enumerate() {
            if i != j && !stop_dupe.contains(&(i, j)) && !stop_dupe.contains(&(j, i)) {
                stop_dupe.insert((i, j));
                stop_dupe.insert((j, i));
                if value > 0.0 {
                    contact_list.push([i as f64, j as f64, value]);
                } else {
                    zero_ind.push(count);
                }
                count += 1;
            }
        }
    }

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
pub fn scale_arr(xyz: &mut [[f64; 3]], min_val: f64, max_val: f64) {
    let flat_min = xyz
        .iter()
        .flat_map(|p| p.iter().copied())
        .fold(f64::INFINITY, f64::min);
    let flat_max = xyz
        .iter()
        .flat_map(|p| p.iter().copied())
        .fold(f64::NEG_INFINITY, f64::max);
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
        writeln!(
            f,
            "{}  {}   {} {}   {}{}{}  {}",
            col1, col2, col3, col4, col5, col6, col7, col8
        )
        .unwrap();
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
pub fn write_log(
    outfile: &str,
    in_file: &str,
    best_alpha: f64,
    rmse: f64,
    best_spearman: f64,
    best_pearson: f64,
) {
    let content = format!(
        "Input file: {}\nConvert factor:: {}\nBest cost  : {}\nBest Spearman correlation Dist vs. Reconstructed Dist  : {}\nBest Pearson correlation Dist vs. Reconstructed Dist  : {}\n",
        in_file, best_alpha, rmse, best_spearman, best_pearson
    );
    fs::write(outfile, content).expect("Failed to write log");
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    const EPS: f64 = 1e-9;

    fn write_fixture(dir: &Path, name: &str, body: &str) -> String {
        let path = dir.join(name);
        fs::write(&path, body).unwrap();
        path.to_str().unwrap().to_string()
    }

    fn sorted(mut values: Vec<f64>) -> Vec<f64> {
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        values
    }

    #[test]
    fn read_matrix_to_list_keeps_upper_triangle_and_records_zero_pairs() {
        let dir = tempdir().unwrap();
        let file = write_fixture(dir.path(), "m.txt", "0 1 2\n1 0 0\n2 0 0\n");

        let (contacts, zero_ind) = read_matrix_to_list(&file);

        assert_eq!(contacts, vec![[0.0, 1.0, 1.0], [0.0, 2.0, 2.0]]);
        assert_eq!(zero_ind, vec![2]);
    }

    #[test]
    fn read_matrix_to_list_drops_all_zero_rows_and_columns() {
        let dir = tempdir().unwrap();
        let file = write_fixture(dir.path(), "m.txt", "0 1 0\n1 0 0\n0 0 0\n");

        let (contacts, zero_ind) = read_matrix_to_list(&file);

        assert_eq!(contacts, vec![[0.0, 1.0, 1.0]]);
        assert!(zero_ind.is_empty());
    }

    #[test]
    fn read_matrix_to_list_returns_nothing_for_an_all_zero_matrix() {
        let dir = tempdir().unwrap();
        let file = write_fixture(dir.path(), "m.txt", "0 0\n0 0\n");

        let (contacts, zero_ind) = read_matrix_to_list(&file);

        assert!(contacts.is_empty());
        assert!(zero_ind.is_empty());
    }

    #[test]
    fn read_matrix_to_list_ignores_blank_lines_and_unparsable_values() {
        let dir = tempdir().unwrap();
        let file = write_fixture(dir.path(), "m.txt", "0 1\n\n1 NA\n");

        let (contacts, zero_ind) = read_matrix_to_list(&file);

        assert_eq!(contacts, vec![[0.0, 1.0, 1.0]]);
        assert!(zero_ind.is_empty());
    }

    #[test]
    fn read_data_remaps_contacts_onto_a_compact_index_range() {
        let dir = tempdir().unwrap();
        let file = write_fixture(dir.path(), "m.txt", "0 1 2\n1 0 0\n2 0 0\n");

        let (contacts, point_map, zero_ind) = read_data(&file);

        assert_eq!(point_map.len(), 3);
        let mut compact: Vec<usize> = point_map.values().copied().collect();
        compact.sort_unstable();
        assert_eq!(compact, vec![0, 1, 2]);

        for c in &contacts {
            assert!((c[0] as usize) < 3);
            assert!((c[1] as usize) < 3);
            assert_ne!(c[0], c[1]);
        }
        assert_eq!(
            sorted(contacts.iter().map(|c| c[2]).collect()),
            vec![1.0, 2.0]
        );
        assert_eq!(zero_ind, vec![2]);
    }

    #[test]
    fn strip_file_collapses_whitespace_into_single_spaces() {
        let dir = tempdir().unwrap();
        let file = write_fixture(dir.path(), "ragged.txt", "1\t\t2   3\n  4 5\t6  \n");

        let out = strip_file(&file);

        assert_eq!(out, format!("{}.stripped", file));
        assert_eq!(fs::read_to_string(&out).unwrap(), "1 2 3\n4 5 6\n");
    }

    #[test]
    fn scale_arr_maps_the_global_extremes_onto_the_requested_bounds() {
        let mut xyz = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];

        scale_arr(&mut xyz, -10.0, 10.0);

        assert_eq!(xyz, vec![[-10.0, -10.0, -10.0], [10.0, 10.0, 10.0]]);
    }

    #[test]
    fn scale_arr_preserves_relative_spacing() {
        let mut xyz = vec![[0.0, 1.0, 2.0], [3.0, 4.0, 4.0]];

        scale_arr(&mut xyz, 0.0, 1.0);

        assert!((xyz[0][0] - 0.0).abs() < EPS);
        assert!((xyz[0][1] - 0.25).abs() < EPS);
        assert!((xyz[0][2] - 0.5).abs() < EPS);
        assert!((xyz[1][0] - 0.75).abs() < EPS);
        assert!((xyz[1][1] - 1.0).abs() < EPS);
    }

    #[test]
    fn write_pdb_emits_one_atom_per_bead_plus_conect_and_end_records() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("chr.pdb");
        let positions = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];

        write_pdb(&positions, out.to_str().unwrap());

        let text = fs::read_to_string(&out).unwrap();
        assert_eq!(text.lines().filter(|l| l.starts_with("ATOM")).count(), 3);
        assert_eq!(text.lines().filter(|l| l.starts_with("CONECT")).count(), 3);
        assert_eq!(text.lines().last().unwrap(), "END");
        assert!(text.contains("   1.000   2.000   3.000"));
    }

    #[test]
    fn write_output_scales_coordinates_into_the_pdb_range() {
        let dir = tempdir().unwrap();
        let base = dir.path().join("chr");
        let positions = [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [4.0, 4.0, 4.0]];

        write_output(base.to_str().unwrap(), &positions);

        let text = fs::read_to_string(format!("{}.pdb", base.to_str().unwrap())).unwrap();
        let coords: Vec<f64> = text
            .lines()
            .filter(|l| l.starts_with("ATOM"))
            .flat_map(|l| {
                let t: Vec<&str> = l.split_whitespace().collect();
                let n = t.len();
                [t[n - 5], t[n - 4], t[n - 3]]
                    .map(|v| v.parse::<f64>().unwrap())
                    .to_vec()
            })
            .collect();

        assert_eq!(coords.len(), 9);
        let lo = coords.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = coords.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!((lo + 10.0).abs() < 1e-3);
        assert!((hi - 10.0).abs() < 1e-3);
    }

    #[test]
    fn write_log_records_the_input_file_and_every_metric() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("run.log");

        write_log(out.to_str().unwrap(), "chr21.txt", 0.3, 1.25, 0.91, 0.87);

        let text = fs::read_to_string(&out).unwrap();
        assert!(text.contains("Input file: chr21.txt"));
        assert!(text.contains("Convert factor:: 0.3"));
        assert!(text.contains("Best cost  : 1.25"));
        assert!(text.contains("Spearman correlation Dist vs. Reconstructed Dist  : 0.91"));
        assert!(text.contains("Pearson correlation Dist vs. Reconstructed Dist  : 0.87"));
    }
}
