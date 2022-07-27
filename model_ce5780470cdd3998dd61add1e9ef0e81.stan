
    data {
        int<lower=0> N;
        array[N] real x;
    }
    parameters {
        real mu;
        real<lower=0> sigma;
    }
    model {
        target += normal_lpdf(mu | 0, 2);  // more informative prior
        target += normal_lpdf(x | mu, sigma);
    }

