class states {

 public:
  double X[NG][RG+1][DS][ART_ST];

  #if fullmod
  double trans[NG][RG][DS][ART_ST];
  double inc[NG][RG];
  double art_init[NG][RG+1][DS];
  #endif

  size_t art_index[DS];
  
  states();

  /* operators */
  states operator+(const states& a);
  states operator*(const double& c);

  states & operator+=(const states& a);

};
