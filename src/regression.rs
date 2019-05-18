use crate::geneobject::{Dataset, QTL};
use rand::Rng;
use rayon::prelude::*;

pub struct RegResult {
    lrs: f64,
    additive: f64,
    dominance: Option<f64>,
}

fn permuted<T>(data: &[T]) -> Vec<T>
where
    T: Copy,
{
    let mut result = data.to_owned();
    let n = data.len();

    for ix in 0..n {
        let j = rand::thread_rng().gen_range(0, n);
        result.swap(ix, j);
    }

    result
}

fn permuted_mut<T>(data: &mut Vec<T>) {
    let n = data.len();
    for ix in 0..n {
        let j = rand::thread_rng().gen_range(0, n);
        data.swap(ix, j);
    }
}

pub fn pvalue(lrs: f64, permutations: &[f64]) -> f64 {
    let n = permutations.len();
    let mut i = 0;

    let mut temp = Vec::from(permutations);
    temp.sort_by(|x, y| x.partial_cmp(y).unwrap());
    for val in temp.iter() {
        if val > &lrs {
            break;
        }
        i += 1;
    }

    if i == n {
        0.0
    } else if i == 0 {
        1.0
    } else {
        1.0 - ((i as f64) / (n as f64))
    }
}

// TODO: add support for variance and control
// TODO: add support for providing a list of strain names to include
pub fn regression(
    dataset: &Dataset,
    traits: &[f64],
    strains: &[String],
    control: Option<&str>,
) -> Vec<QTL> {
    //
    let mut result = Vec::with_capacity(dataset.n_loci());

    let strain_ixs = dataset.strain_indices(strains);

    let control_geno: Option<Vec<_>> = control.map(|c| {
        dataset
            .genome
            .find_locus(c)
            .unwrap()
            .genotypes_subset(&strain_ixs)
    });

    if control != None && control_geno == None {
        panic!("Control could not be found in loci list");
    }

    for loci in dataset.genome.iter() {
        for locus in loci.iter() {
            let genotypes = locus.genotypes_subset(&strain_ixs);

            let reg_result = match &control_geno {
                None => {
                    if dataset.dominance {
                        let dominance = locus.dominance_subset(&strain_ixs);
                        regression_3n(traits, &genotypes, &dominance, false)
                    } else {
                        regression_2n(traits, &genotypes)
                    }
                }
                Some(c) => {
                    if dataset.dominance {
                        panic!("reaper: no composite regression for intercross");
                    } else {
                        regression_3n(traits, &genotypes, &c, false)
                    }
                }
            };

            result.push(QTL {
                lrs: reg_result.lrs,
                additive: reg_result.additive,
                dominance: reg_result.dominance,
                marker: locus.marker.clone(),
            })
        }
    }

    result
}

pub fn permutation(
    dataset: &Dataset,
    traits: &[f64],
    strains: &[String],
    n_perms: usize,
    threads: usize,
) -> Vec<f64> {
    let threads = std::cmp::max(threads, 1);
    // let lrs_thresh = -1.0;
    // let top_n = 10;

    let strain_ixs = dataset.strain_indices(strains);

    let mut vecs = Vec::with_capacity(threads);
    vecs.par_extend((0..threads).into_par_iter().map(|_| {
        let mut temp_vec = Vec::with_capacity(n_perms / 4);
        let mut p_traits = permuted(traits);
        (0..(n_perms / threads)).for_each(|_| {
            let mut lrs_max = 0.0;
            let mut genotypes = vec![0.0; strain_ixs.len()];

            for loci in dataset.genome.iter() {
                for locus in loci.iter() {
                    locus.genotypes_subindices(&strain_ixs, &mut genotypes);
                    let reg_result = regression_2n(&p_traits, &genotypes);
                    lrs_max = reg_result.lrs.max(lrs_max);
                }
            }
            temp_vec.push(lrs_max);

            permuted_mut(&mut p_traits);
        });
        temp_vec.into_iter()
    }));
    let mut lrs_vec: Vec<_> = vecs.into_iter().flatten().collect();

    lrs_vec.sort_by(|x, y| x.partial_cmp(y).unwrap());
    lrs_vec
}
// `traits` corresponds to `YY`
// `genotypes` corresponds to `XX`
fn regression_2n(traits: &[f64], genotypes: &[f64]) -> RegResult {
    let mut sig_y = 0.0;
    let mut sig_yy = 0.0;

    let mut sig_x = 0.0;
    let mut sig_xx = 0.0;
    let mut sig_xy = 0.0;

    let n_strains = traits.len();
    let n = n_strains as f64;

    for ix in 0..traits.len() {
        let temp_trait = traits[ix];
        let temp_geno = genotypes[ix];

        sig_y += temp_trait;
        sig_yy += temp_trait * temp_trait;
        sig_xy += temp_trait * temp_geno;

        sig_x += temp_geno;
        sig_xx += temp_geno * temp_geno;
    }

    let d = sig_xx - sig_x * sig_x / n;
    let tss = sig_yy - (sig_y * sig_y) / n;

    let a = (sig_xx * sig_y - sig_x * sig_xy) / (n * d);
    let mut b = (sig_xy - (sig_x * sig_y / n)) / d;

    let rss =
        sig_yy + a * (n * a - 2.0 * sig_y) + b * (2.0 * a * sig_x + b * sig_xx - 2.0 * sig_xy);

    let mut lrs = n * (tss / rss).ln();

    if lrs.is_nan() || lrs < 0.0 {
        b = 0.0;
        lrs = 0.0;
    }

    RegResult {
        lrs,
        additive: b,
        dominance: None,
    }
}

fn regression_2n_variance(traits: &[f64], genotypes: &[f64], variance: &[f64]) -> RegResult {
    let mut sigYV = 0.0;
    let mut sigYYV = 0.0;
    let mut sigXV = 0.0;
    let mut sigXXV = 0.0;
    let mut sigXYV = 0.0;

    let mut sig1V = 0.0;

    let n_strains = traits.len();

    for ix in 0..traits.len() {
        let temp0 = 1.0 / variance[ix];
        let temp1 = traits[ix];
        let temp2 = genotypes[ix];
        sig1V += temp0;
        let temp = temp1 * temp0;
        sigYV += temp;
        sigYYV += temp1 * temp;
        sigXYV += temp * temp2;
        let temp = temp2 * temp0;
        sigXV += temp;
        sigXXV += temp * temp2;
    }

    let d = sigXXV - sigXV * sigXV / sig1V;
    let tss = sigYYV - (sigYV * sigYV) / sig1V;
    let a = (sigXXV * sigYV - sigXV * sigXYV) / (sig1V * d);
    let mut b = (sigXYV - (sigXV * sigYV / sig1V)) / d;
    let rss =
        sigYYV + a * (sig1V * a - 2.0 * sigYV) + b * (2.0 * a * sigXV + b * sigXXV - 2.0 * sigXYV);
    let mut lrs = (n_strains as f64) * (tss / rss).ln();

    if lrs.is_nan() || lrs < 0.0 {
        b = 0.0;
        lrs = 0.0;
    }

    RegResult {
        lrs,
        additive: b,
        dominance: None,
    }
}

fn regression_3n(traits: &[f64], genotypes: &[f64], controls: &[f64], diff: bool) -> RegResult {
    let mut sigC = 0.0;
    let mut sigX = 0.0;
    let mut sigY = 0.0;
    let mut sigCC = 0.0;
    let mut sigXX = 0.0;
    let mut sigYY = 0.0;
    let mut sigXC = 0.0;
    let mut sigCY = 0.0;
    let mut sigXY = 0.0;

    let n_strains = traits.len();
    let n = n_strains as f64;

    for ix in 0..traits.len() {
        let a = controls[ix];
        let b = genotypes[ix];
        let y = traits[ix];
        sigC += a;
        sigX += b;
        sigY += y;
        sigCC += a * a;
        sigXX += b * b;
        sigYY += y * y;
        sigXC += a * b;
        sigCY += y * a;
        sigXY += y * b;
    }

    let temp0 = sigXC * sigXC - sigCC * sigXX;
    let temp1 = sigC * sigXX - sigX * sigXC;
    let temp2 = sigX * sigCC - sigC * sigXC;
    let temp3 = sigX * sigX - n * sigXX;
    let temp4 = n * sigXC - sigC * sigX;
    let temp5 = sigC * sigC - n * sigCC;
    let temp6 = temp4 * sigXC + temp2 * sigX + temp5 * sigXX;

    let betak = (temp0 * sigY + temp1 * sigCY + temp2 * sigXY) / temp6;
    let mut betac = (temp1 * sigY + temp3 * sigCY + temp4 * sigXY) / temp6;
    let mut betax = (temp2 * sigY + temp4 * sigCY + temp5 * sigXY) / temp6;

    let ssf = sigYY
        + betac * (betac * sigCC - 2.0 * sigCY)
        + betax * (betax * sigXX - 2.0 * sigXY)
        + 2.0 * betac * betax * sigXC
        + betak * (n * betak + 2.0 * betac * sigC + 2.0 * betax * sigX - 2.0 * sigY);

    let ssr = if diff {
        let d = sigCC - sigC * sigC / n;
        let a = (sigCC * sigY - sigC * sigCY) / (n * d);
        let b = (sigCY - (sigC * sigY / n)) / d;
        sigYY + a * (n * a - 2.0 * sigY) + b * (2.0 * a * sigC + b * sigCC - 2.0 * sigCY)
    } else {
        sigYY - (sigY * sigY) / n
    };

    let mut lrs = n * (ssr / ssf).ln();
    if lrs.is_nan() || lrs < 0.0 {
        betax = 0.0;
        lrs = 0.0;
        // NOTE: in the old implementation it is `betak`, not `betac`, that is set to 0.0 here, but `betak` is not used later, so I assume it's a mistake!
        betac = 0.0;
    }

    RegResult {
        lrs,
        additive: betax,
        dominance: Some(betac),
    }
}

// this one will require a bit more work since it actually uses matrices!
/*
fn regression_3n(
    traits: &[f64],
    genotypes: &[f64],
    controls: &[f64],
    variance: &[f64],
    diff: bool,
) -> RegResult {


    let mut sig1V = 0.0;
    let mut sigYV = 0.0;
    let mut sigXV =0.0;
    let mut sigCV =0.0;
    let mut sigXXV = 0.0;
    let mut sigYYV = 0.0;
    let mut sigCCV = 0.0;
    let mut sigXYV = 0.0;
    let mut sigXCV = 0.0;
    let mut sigCYV = 0.0;


    let n_strains = traits.len();
    let n = n_strains as f64;

    for ix in 0..traits.len() {
        let c = controls[ix];
        let x = genotypes[ix];
        let y = traits[ix];
        let v = 1.0/variance[ix];
        sig1V += v;
        sigYV += y*v;
        sigXV += x*v;
        sigCV += c*v;
        sigXXV += x*x*v;
        sigYYV += y*y*v;
        sigCCV += c*c*v;
        sigXYV += x*y*v;
        sigXCV += x*c*v;
        sigCYV += c*y*v;
    }

    aa = square(3);
    aa[0][0] = sig1V;
    aa[1][1] = sigXXV;
    aa[2][2] = sigCCV;
    aa[0][1] = aa[1][0] = sigXV;
    aa[0][2] = aa[2][0] = sigCV;
    aa[1][2] = aa[2][1] = sigXCV;

    inverse(aa,3);

    betak = aa[0][0]*sigYV + aa[0][1]*sigXYV + aa[0][2]*sigCYV;
    betax = aa[1][0]*sigYV + aa[1][1]*sigXYV + aa[1][2]*sigCYV;
    betac = aa[2][0]*sigYV + aa[2][1]*sigXYV + aa[2][2]*sigCYV;
    ssf = sigYYV+betax*(betax*sigXXV-2*sigXYV)+betac*(betac*sigCCV-2*sigCYV)+ 2*betax*betac*sigXCV+betak*(sig1V*betak+2*betax*sigXV+2*betac*sigCV-2*sigYV);
    if (diff != 0){
        D = sigCCV - sigCV*sigCV/sig1V;
        a = (sigCCV*sigYV - sigCV*sigCYV)/(sig1V*D);
        b = (sigCYV - (sigCV*sigYV/sig1V))/D;
        ssr = sigYYV + a*(sig1V*a-2.0*sigYV) + b*(2.0*a*sigCV+b*sigCCV-2.0*sigCYV);
    }
    else
        ssr = sigYYV - (sigYV*sigYV)/sig1V;

    LRS = n*log(ssr/ssf);
    if (isnan(LRS) || (LRS < 0)){
        betax = 0.0;
        betac = 0.0;
        LRS = 0.0;
    }
    result[0] = LRS;
    result[1] = betax;
    result[2] = betac;
    freesquare(aa,3);
    return 1;

    RegResult {
        lrs,
        additive: betax,
        dominance: Some(betac),
    }
}

*/
