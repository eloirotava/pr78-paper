use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

const R: f64 = 8.314_462_618; // J/mol/K

// CO2 (valores típicos; se quiser, depois você substitui por dados gerados do CoolProp)
const TC: f64 = 304.1282;     // K
const PC: f64 = 7.3773e6;     // Pa
const OMEGA: f64 = 0.22394;   // -

fn kappa(omega: f64) -> f64 {
    0.37464 + 1.54226 * omega - 0.26992 * omega * omega
}

fn alpha(t: f64, tc: f64, omega: f64) -> f64 {
    let tr = t / tc;
    let k = kappa(omega);
    let s = 1.0 + k * (1.0 - tr.sqrt());
    s * s
}

fn a_c(tc: f64, pc: f64) -> f64 {
    0.45724 * R * R * tc * tc / pc
}

fn b_i(tc: f64, pc: f64) -> f64 {
    0.07780 * R * tc / pc
}

fn pr_pressure_z_co2(rho: f64, t: f64) -> Result<(f64, f64), String> {
    if rho <= 0.0 || t <= 0.0 {
        return Err("rho e T precisam ser positivos".into());
    }

    let a = a_c(TC, PC) * alpha(t, TC, OMEGA);
    let b = b_i(TC, PC);

    let v = 1.0 / rho;
    let vb = v - b;
    if vb <= 0.0 {
        return Err("estado inválido: v-b<=0 (rho alto demais)".into());
    }

    // denom = v(v+b)+b(v-b) = v^2 + 2bv - b^2
    let denom = v * v + 2.0 * b * v - b * b;
    if denom <= 0.0 {
        return Err("estado inválido: denom<=0".into());
    }

    let p = R * t / vb - a / denom;
    let z = p / (rho * R * t);
    Ok((p, z))
}

fn main() -> Result<(), String> {
    std::fs::create_dir_all("out").map_err(|e| e.to_string())?;

    // Grade perto do crítico
    let t_min = TC + 0.0;
    let t_max = TC + 10.0;
    let n_t: usize = 81;

    // faixa de densidade molar (CO2 crítico ~ 10.6 kmol/m3 = 10600 mol/m3)
    let rho_min = 2000.0;
    let rho_max = 16000.0;
    let n_rho: usize = 120;

    let f = File::create("out/pr78_co2_grid.csv").map_err(|e| e.to_string())?;
    let mut w = BufWriter::new(f);

    writeln!(w, "T_K,rho_mol_m3,P_Pa,Z_PR").map_err(|e| e.to_string())?;

    for it in 0..n_t {
        let t = t_min + (t_max - t_min) * (it as f64) / ((n_t - 1) as f64);
        for ir in 0..n_rho {
            let rho = rho_min + (rho_max - rho_min) * (ir as f64) / ((n_rho - 1) as f64);

            match pr_pressure_z_co2(rho, t) {
                Ok((p, z)) => {
                    writeln!(w, "{:.6},{:.6},{:.12e},{:.12}", t, rho, p, z)
                        .map_err(|e| e.to_string())?;
                }
                Err(_) => {
                    // Se der estado inválido, você pode pular ou escrever NaN:
                    writeln!(w, "{:.6},{:.6},NaN,NaN", t, rho).map_err(|e| e.to_string())?;
                }
            }
        }
    }

    println!("Gerado: out/pr78_co2_grid.csv");
     
    let mut acc = 0.0_f64; // evita otimização boba
    let t0 = Instant::now();

    for it in 0..n_t {
        let t = t_min + (t_max - t_min) * (it as f64) / ((n_t - 1) as f64);
        for ir in 0..n_rho {
            let rho = rho_min + (rho_max - rho_min) * (ir as f64) / ((n_rho - 1) as f64);
            if let Ok((_p, z)) = pr_pressure_z_co2(rho, t) {
                acc += z;
            }
        }
    }

    let dt = t0.elapsed();
    let n = (n_t * n_rho) as f64;
    println!("Rust PR78: pontos = {}, tempo = {:?}, pontos/s = {:.3e}, checksum = {:.6}",
             n_t*n_rho, dt, n / dt.as_secs_f64(), acc);

    Ok(())
    
}
