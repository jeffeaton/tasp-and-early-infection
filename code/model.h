class states;
states rk4(states y, double dt, const struct parameters *);
states euler(states y, const double dt, const struct parameters *);
states grad(states y, const struct parameters *);

// Parameters //

extern const double alpha;
extern const double nu;
extern const double mu;
extern const double theta_g;
extern const double sigma[NG][DS][ART_ST-8];
extern const double xi[NG][DS];
extern const double eta[NG][DS][ART_ST-9];
extern const double eta_reinit[NG][DS][ART_ST-9];
extern const double lambda_reinit[NG][DS];
extern const double lambda_reinit_late[NG][DS];
