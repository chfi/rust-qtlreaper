pub struct RegResult {
    lrs: f64,
    additive: f64,
    dominance: Option<f64>,
}

// `traits` corresponds to `YY`
// `genotypes` corresponds to `XX`
pub fn regression_2n(traits: Vec<f64>, genotypes: Vec<f64>) -> RegResult {
    let mut sig_y = 0.0;
    let mut sig_yy = 0.0;

    let mut sig_x = 0.0;
    let mut sig_xx = 0.0;
    let mut sig_xy = 0.0;

    let n_strains = traits.len();

    for ix in 0..traits.len() {
        let temp_trait = traits[ix];
        let temp_geno = genotypes[ix];

        sig_y += temp_trait;
        sig_yy += temp_trait * temp_trait;
        sig_yy += temp_trait * temp_geno;

        sig_x += temp_geno;
        sig_xx += temp_geno * temp_geno;
    }

    let n = n_strains as f64;

    let d = sig_xx - sig_x * sig_x / n;
    let tss = sig_yy - (sig_y * sig_y) / n;

    let a = (sig_xx * sig_y - sig_x * sig_xy) / (n * d);
    let mut b = (sig_xy - (sig_x * sig_y / n)) / d;

    let rss =
        sig_yy + a * (n * a - 2.0 * sig_y) + b * (2.0 * a * sig_x + b * sig_xx - 2.0 * sig_xy);

    let mut lrs = n * (tss / rss).ln();

    if lrs == std::f64::NAN || lrs < 0.0 {
        b = 0.0;
        lrs = 0.0;
    }

    RegResult {
        lrs,
        additive: b,
        dominance: None,
    }
}

pub fn regression_2n_variance(
    traits: Vec<f64>,
    genotypes: Vec<f64>,
    variance: Vec<f64>,
) -> RegResult {
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

    if lrs == std::f64::NAN || lrs < 0.0 {
        b = 0.0;
        lrs = 0.0;
    }

    RegResult {
        lrs,
        additive: b,
        dominance: None,
    }
}

pub fn regression_3n(
    traits: Vec<f64>,
    genotypes: Vec<f64>,
    controls: Vec<f64>,
    diff: bool,
) -> RegResult {
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

    let mut betak = (temp0 * sigY + temp1 * sigCY + temp2 * sigXY) / temp6;
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
    if lrs == std::f64::NAN || lrs < 0.0 {
        betax = 0.0;
        betak = 0.0;
        lrs = 0.0;
    }

    RegResult {
        lrs,
        additive: betax,
        dominance: Some(betac),
    }
}

// this one will require a bit more work since it actually uses matrices!
/*
pub fn regression_3n(
    traits: Vec<f64>,
    genotypes: Vec<f64>,
    controls: Vec<f64>,
    variance: Vec<f64>,
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
