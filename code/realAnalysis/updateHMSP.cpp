// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;
using namespace Eigen;

// [[Rcpp::export]]
void check_matrix_size(const NumericMatrix& A) {
  // 获取矩阵的行数和列数
  int n_rows = A.nrow();
  int n_cols = A.ncol();
  
  // 打印矩阵的维度
  Rcpp::Rcout << "Matrix A: " << n_rows << " x " << n_cols << std::endl;
  
  // 计算并打印矩阵的大小 (单位: MB)
  double size_in_mb = (n_rows * n_cols * sizeof(double)) / (1024.0 * 1024.0);
  Rcpp::Rcout << "Size of Matrix A (in MB): " << size_in_mb << " MB" << std::endl;
}


// [[Rcpp::export]]
NumericMatrix f_madd(NumericMatrix tmm, NumericMatrix tm22){
  
  const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
  const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));
  
  Eigen::MatrixXd sum = ttm+ttm2;
  return(wrap(sum));
}

// [[Rcpp::export]]
NumericMatrix f_mmult(NumericMatrix tmm, NumericMatrix tm22){
  
  const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
  const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));
  
  Eigen::MatrixXd prod = ttm*ttm2;
  return(wrap(prod));
}

// [[Rcpp::export]]
NumericMatrix f_mvnorm(int n, NumericVector mu, NumericMatrix sigma) {
  arma::mat sigma_temp= Rcpp::as<arma::mat>(sigma);
  arma::vec mu_temp = Rcpp::as<arma::vec>(mu);
  int ncols = sigma_temp.n_cols;
  mat Y = randn(n, ncols);
  return(wrap((repmat(mu_temp, 1, n).t()+Y * chol(sigma_temp)).t()));
}


// [[Rcpp::export]]
NumericVector f_invgauss(int n, double mu, double lambda){
  
  NumericVector random_vector(n);
  double z,y,x,u;
  
  for(int i=0; i<n; ++i){
    z=rnorm(1,0,1)[0];
    y=z*z;
    x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
    u=runif(1,0,1)[0];
    if(u <= mu/(mu+x)){
      random_vector(i)=x;
    }else{
      random_vector(i)=mu*mu/x;
    };
  }
  return(random_vector);
}


// [[Rcpp::export]]
NumericMatrix f_inv(NumericMatrix tmm){
  
  const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
  
  Eigen::MatrixXd inv = ttm.inverse();
  return(wrap(inv));
}


// [[Rcpp::export]]
NumericVector rdirichlet(int n, NumericVector alpha) {
  int k = alpha.size();          // Dirichlet 分布的维度
  NumericVector samples(k);      // 存储伽马分布随机数
  double sum_samples = 0.0;      // 用于归一化
  
  // 生成每个维度的伽马分布随机数
  for (int i = 0; i < k; ++i) {
    if (alpha[i] <= 0) {
      stop("All alpha values must be positive.");
    }
    samples[i] = R::rgamma(alpha[i], 1.0);  // 调用 R 的伽马分布生成器
    sum_samples += samples[i];
  }
  
  // 归一化使和为 1
  for (int i = 0; i < k; ++i) {
    samples[i] /= sum_samples;
  }
  
  return samples;
}


// [[Rcpp::export]]
double f_det(const NumericMatrix A) {
  // 将 NumericMatrix 转换为 Eigen::Map<Eigen::MatrixXd>
    const Eigen::Map<Eigen::MatrixXd> mat(as<Eigen::Map<Eigen::MatrixXd>>(A));
  
  // 检查是否为方阵
  if (mat.rows() != mat.cols()) {
    stop("The matrix must be square to compute the determinant.");
  }
  
  // 使用 Eigen 的 determinant() 计算行列式
  return mat.determinant();
}


// [[Rcpp::export]]
NumericMatrix f_transpose(const NumericMatrix A) {
  int nrow = A.nrow();
  int ncol = A.ncol();
  
  // 创建一个新的矩阵，大小为转置后的尺寸
  NumericMatrix At(ncol, nrow);
  
  // 手动交换行列元素
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      At(j, i) = A(i, j);  // 交换行列位置
    }
  }
  
  return At;
}


// [[Rcpp::export]]
double f_maxValue(const NumericMatrix& A) {
  double max_val = A(0, 0);  // 初始化最大值为第一个元素
  int nrow = A.nrow();
  int ncol = A.ncol();
  
  // 遍历矩阵所有元素，找到最大值
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (A(i, j) > max_val) {
        max_val = A(i, j);
      }
    }
  }
  
  return max_val;
}

// [[Rcpp::export]]
Rcpp::NumericVector f_order(arma::vec x) {
  return(Rcpp::as<Rcpp::NumericVector>(wrap(arma::sort_index( x )+1)) );
}


// [[Rcpp::export]]
NumericVector f_sample( NumericVector x, 
                        int size,
                        bool replace, 
                        NumericVector prob = NumericVector::create()
){
  NumericVector ret = Rcpp::sample(x, size, replace, prob) ;
  return ret ;
};


// [[Rcpp::export]]
NumericVector f_likeli(NumericVector Bbeta, NumericVector y, NumericVector sigma, int N, int NT, int NS) {
  NumericVector L_out(NS*N*NT);
  for(int i=0; i<NS; ++i){
    for(int j=0; j<(N*NT); ++j){
      L_out[i*N*NT+j]=exp(-0.5*(log(sigma[i])+pow((y[j]-Bbeta[i*N*NT+j]),2)/sigma[i]));
      // Rcpp::Rcout << "The likelihood value is: " << L_out[i*N*NT+j] << std::endl;
    }
    // Rcpp::Rcout << "Likelihood Elements: ";
    // for (int j=0; j<(N*NT); ++j) {
      //   Rcpp::Rcout << L_out[i*N*NT+j] << " ";
      // }
    // Rcpp::Rcout << std::endl;  // 换行
  }
  return L_out;
}


// [[Rcpp::export]]
List f_p(NumericVector tran, NumericVector tau, int N, int NT, int NS){
  NumericVector p_out(N*NT*NS*NS);
  NumericVector p0_out(NS);
  for(int i=1; i<NT; ++i){
    for(int j=0; j<N; ++j){
      for(int k=0; k<NS; ++k){
        for(int l=0; l<NS; ++l){
          p_out[i*N*NS*NS+j*NS*NS+k*NS+l]=1;
          if(l==0){
            p_out[i*N*NS*NS+j*NS*NS+k*NS+l]=exp(tran[i*N*NS*NS+j*NS*NS+k*NS+l])/(1+exp(tran[i*N*NS*NS+j*NS*NS+k*NS+l]));
          }else if(l==(NS-1)){
            for(int m=0; m<(NS-1); ++m){
              p_out[i*N*NS*NS+j*NS*NS+k*NS+l] *= 1/(1+exp(tran[i*N*NS*NS+j*NS*NS+k*NS+m]));
            }
          }else{
            for(int m=0; m<l; ++m){
              p_out[i*N*NS*NS+j*NS*NS+k*NS+l] *= 1/(1+exp(tran[i*N*NS*NS+j*NS*NS+k*NS+m]));
            }
            p_out[i*N*NS*NS+j*NS*NS+k*NS+l] *= exp(tran[i*N*NS*NS+j*NS*NS+k*NS+l])/(1+exp(tran[i*N*NS*NS+j*NS*NS+k*NS+l]));
          }
        }
      }
    }
  }
  for(int i=0; i<NS; ++i){
    p0_out[i]=1;
    if(i==0){
      p0_out[i]=exp(tau[i])/(1+exp(tau[i]));
    }else if(i==(NS-1)){
      for(int m=0; m<(NS-1); ++m){
        p0_out[i] *= 1/(1+exp(tau[m]));
      }
    }else{
      for(int m=0; m<i; ++m){
        p0_out[i] *= 1/(1+exp(tau[m]));
      }
      p0_out[i] *= exp(tau[i])/(1+exp(tau[i]));
    }
  }
  return List::create(Rcpp::Named("p") = p_out,
                      Rcpp::Named("p0") = p0_out);
}


// [[Rcpp::export]]
NumericVector f_tran(NumericVector d, NumericVector zeta, NumericVector phi, int N, int NT, int NS, int ND){
  NumericVector tran_out(N*NT*NS*NS);
  for(int i=0; i<NT; ++i){
    for(int j=0; j<N; ++j){
      NumericVector d_temp(ND);
      for(int q=0; q<ND; q++){
        d_temp[q]=d[q*N*NT+i*N+j];
      }
      for(int k=0; k<NS; ++k){
        for(int l=0; l<NS; ++l){
          tran_out[i*N*NS*NS+j*NS*NS+k*NS+l]=zeta[k*NS+l]+sum(d_temp*phi);
        }
      }
    }
  }
  return tran_out;
}



// [[Rcpp::export]]
NumericVector f_s(NumericVector L, NumericVector p, NumericVector p0, int N, int NT, int NS){
  NumericVector Q_forward(N*NT*NS);
  NumericVector Q_back(N*NT*NS);
  NumericVector Q(N*NT*NS);
  NumericVector s_out(N*NT);
  for(int i=0; i<N; ++i){
    for(int j=0; j<NS; ++j){
      Q_forward[i*NS+j]=L[j*N*NT+i]*p0[j];
      Q_back[(NT-1)*N*NS+i*NS+j]=1;
    }
    for(int k=1; k<NT; ++k){
      for(int j=0; j<NS; ++j){
        for(int l=0; l<NS; ++l){
          Q_forward[k*N*NS+i*NS+j] += Q_forward[(k-1)*N*NS+i*NS+l]*p[k*N*NS*NS+i*NS*NS+l*NS+j]*L[j*N*NT+k*N+i];
        }
      }
    }
    for(int k=(NT-2); k>(-1); --k){
      for(int j=0; j<NS; ++j){
        for(int l=0; l<NS; ++l){
          Q_back[k*N*NS+i*NS+j] += Q_back[(k+1)*N*NS+i*NS+l]*p[(k+1)*N*NS*NS+i*NS*NS+j*NS+l]*L[l*N*NT+(k+1)*N+i];
        }
      }
    }
    NumericVector sample(NS);
    for(int j=0; j<NS; ++j){
      sample[j]=j;
    }
    for(int k=0; k<NT; ++k){
      NumericVector prob_temp(NS);
      for(int j=0; j<NS; ++j){
        Q[k*N*NS+i*NS+j]=Q_forward[k*N*NS+i*NS+j]*Q_back[k*N*NS+i*NS+j];
        prob_temp[j]=Q[k*N*NS+i*NS+j];
      }
      
      // Rcpp::Rcout << "Sample Prob Elements: ";
      // for (int j=0; j<(NS); ++j) {
        //   Rcpp::Rcout << prob_temp[j] << " ";
        // }
      // Rcpp::Rcout << std::endl;  // 换行
      // Rcpp::Rcout << "The sample prob value is: " << prob_temp << std::endl;
      
      s_out[k*N+i]=f_sample(sample,1,FALSE,prob_temp)[0];
      // Rcpp::Rcout << "The sample result is: " << s_out[k*N+i] + 1 << std::endl;
    }
  }
  return s_out;
}


// [[Rcpp::export]]
NumericVector f_trprob(NumericVector s, int N, int NT, int NS, NumericVector ele){
  
  NumericVector p(NS*NS);
  NumericVector tr_times(NS*NS);
  // Rcpp::Rcout << "The tr times vector is: " << tr_times << std::endl;
  for(int t=0; t<(N*NT - 1); ++t){
    for(int i=0; i<NS; ++i){
      for(int j=0; j<NS; ++j){
        if (s(t) == i && s(t+1) == j){
          tr_times(j*NS+i) += 1;
        }
      }
    }
  }
  // Rcpp::Rcout << "The tr times vector is: " << tr_times << std::endl;
  NumericVector vec_dir(NS);
  NumericVector alpha(NS);
  for(int i=0; i<NS; ++i){
    for(int j=0; j<NS; ++j){
      alpha(j) = tr_times(j*NS+i) + 1;
      if (i == j){
        alpha(j) += ele(j);
      }
    }
    vec_dir = rdirichlet(1, alpha);
    for(int j=0; j<NS; ++j){
      p(j*NS+i) = vec_dir(j);
    }
  }
  return p;
}


// [[Rcpp::export]]
NumericVector f_s_hmm(NumericVector L, NumericVector p, NumericVector p0, int N, int NT, int NS){
  NumericVector Q_forward(N*NT*NS);
  NumericVector Q_back(N*NT*NS);
  NumericVector Q(N*NT*NS);
  NumericVector s_out(N*NT);
  for(int i=0; i<N; ++i){
    for(int j=0; j<NS; ++j){
      Q_forward[i*NS+j]=L[j*N*NT+i]*p0[j];
      Q_back[(NT-1)*N*NS+i*NS+j]=1;
    }
    for(int k=1; k<NT; ++k){
      for(int j=0; j<NS; ++j){
        for(int l=0; l<NS; ++l){
          Q_forward[k*N*NS+i*NS+j] += Q_forward[(k-1)*N*NS+i*NS+l]*p[l*NS+j]*L[j*N*NT+k*N+i];
          // Q_forward[k*N*NS+i*NS+j] += Q_forward[(k-1)*N*NS+i*NS+l]*p[j*NS+l]*L[j*N*NT+k*N+i];
        }
      }
    }
    for(int k=(NT-2); k>(-1); --k){
      for(int j=0; j<NS; ++j){
        for(int l=0; l<NS; ++l){
          Q_back[k*N*NS+i*NS+j] += Q_back[(k+1)*N*NS+i*NS+l]*p[j*NS+l]*L[l*N*NT+(k+1)*N+i];
          // Q_back[k*N*NS+i*NS+j] += Q_back[(k+1)*N*NS+i*NS+l]*p[l*NS+j]*L[l*N*NT+(k+1)*N+i];
        }
      }
    }
    NumericVector sample(NS);
    for(int j=0; j<NS; ++j){
      sample[j]=j;
    }
    for(int k=0; k<NT; ++k){
      NumericVector prob_temp(NS);
      for(int j=0; j<NS; ++j){
        Q[k*N*NS+i*NS+j]=Q_forward[k*N*NS+i*NS+j]*Q_back[k*N*NS+i*NS+j];
        prob_temp[j]=Q[k*N*NS+i*NS+j];
      }
      s_out[k*N+i]=f_sample(sample,1,FALSE,prob_temp)[0];
    }
  }
  return s_out;
}


// [[Rcpp::export]]
NumericVector f_tau(NumericVector tau, NumericVector s, NumericVector at_tau, double c_tau, int N, int NS){
  NumericVector tau_out(NS-1);
  NumericVector s1(N);
  for(int i=0; i<N; ++i){
    s1[i]=s[i];
  }
  for(int i=0; i<(NS-1); ++i){
    double sigma_tau = 0;
    double miu_tau = tau[i];
    for(int j=i; j<NS;++j){
      sigma_tau += 0.25*table(s1)[j];
    }
    sigma_tau = 1/(sigma_tau+1);
    double tau_star = rnorm(1,miu_tau,sqrt(c_tau*sigma_tau))[0];
    double ratio_tau = 0;
    for(int j=i; j<NS; ++j){
      ratio_tau += table(s1)[j]*(log(1+exp(miu_tau))-log(1+exp(tau_star)));
    }
    ratio_tau = exp(table(s1)[i]*(tau_star-miu_tau)+ratio_tau-0.5*pow(tau_star,2)+0.5*pow(miu_tau,2));
    double randnum = runif(1)[0];
    if(randnum<ratio_tau){
      at_tau[i] += 1;
      tau_out[i]=tau_star;
    }else{
      tau_out[i]=miu_tau;
    }
  }
  return tau_out;
}


// [[Rcpp::export]]
NumericVector f_phi(NumericVector d, NumericVector zeta, NumericVector phi, NumericVector s, int at_phi, double c_phi, int N, int NT, int ND, int NS){
  NumericVector phi_out(ND);
  NumericMatrix temp_phi(ND,ND);
  NumericMatrix sigma_phi(ND,ND);
  NumericVector miu_phi(ND);
  NumericMatrix phi_star(ND,1);
  for(int i=0; i<NS; ++i){
    for(int j=0; j<(NS-1); ++j){
      for(int k=0; k<(j+1); ++k){
        for(int l=1; l<NT; ++l){
          
          for(int m=0; m<N; ++m){
            if((s[l*N+m]==j)&&(s[(l-1)*N+m]==i)){
              NumericMatrix d_temp1(ND,1);
              NumericMatrix d_temp2(1,ND);
              for(int n=0; n<ND; ++n){
                d_temp1(n,0)=d[n*N*NT+l*N+m];
                d_temp2(0,n)=d[n*N*NT+l*N+m];
              }
              temp_phi = f_madd(temp_phi,f_mmult(d_temp1,d_temp2)*exp(zeta[i*NS+k])/pow((1+exp(zeta[i*NS+k])),2));
            }
          }
        }
      }
    }
    for(int k=0; k<(NS-1); ++k){
      for(int l=1; l<NT; ++l){
        for(int m=0; m<N; ++m){
          if((s[l*N+m]==(NS-1))&&(s[(l-1)*N+m]==i)){
            NumericMatrix d_temp1(ND,1);
            NumericMatrix d_temp2(1,ND);
            for(int n=0; n<ND; ++n){
              d_temp1(n,0)=d[n*N*NT+l*N+m];
              d_temp2(0,n)=d[n*N*NT+l*N+m];
            }
            temp_phi = f_madd(temp_phi,f_mmult(d_temp1,d_temp2)*exp(zeta[i*NS+k])/pow((1+exp(zeta[i*NS+k])),2));
          }
        }
      }
    }
  }
  NumericMatrix one(ND,ND);
  for(int i=0; i<ND; ++i){
    one(i,i)=1;
  }
  sigma_phi = f_inv(f_madd(temp_phi,one));
  miu_phi = phi;
  phi_star = f_mvnorm(1,miu_phi,c_phi*sigma_phi);
  
  double ratio_phi = 0;
  for(int i=0; i<NS; ++i){
    for(int j=0; j<(NS-1); ++j){
      for(int l=1; l<NT; ++l){
        for(int m=0; m<N; ++m){
          if((s[l*N+m]==j)&&(s[(l-1)*N+m]==i)){
            double phixstar = 0;
            double phix = 0;
            for(int n=0; n<ND; ++n){
              phixstar += d[n*N*NT+l*N+m]*phi_star(n,0);
              phix += d[n*N*NT+l*N+m]*miu_phi[n];
            }
            ratio_phi += phixstar-phix;
            for(int k=0; k<(j+1); ++k){
              ratio_phi += log(1+exp(zeta[i*NS+k]+phix))-log(1+exp(zeta[i*NS+k]+phixstar));
            }
          }
        }
      }
    }
    for(int l=1; l<NT; ++l){
      for(int m=0; m<N; ++m){
        
        if((s[l*N+m]==(NS-1))&&(s[(l-1)*N+m]==i)){
          double phixstar = 0;
          double phix = 0;
          for(int n=0; n<ND; ++n){
            phixstar += d[n*N*NT+l*N+m]*phi_star(n,0);
            phix += d[n*N*NT+l*N+m]*miu_phi[n];
          }
          for(int k=0; k<(NS-1); ++k){
            ratio_phi += log(1+exp(zeta[i*NS+k]+phix))-log(1+exp(zeta[i*NS+k]+phixstar));
          }
        }
      }
    }
  }
  double phi2 = 0;
  double phi2star =0;
  for(int i=0; i<ND; ++i){
    phi2 += pow(miu_phi[i],2);
    phi2star += pow(phi_star(i,0),2);
  }
  ratio_phi = exp(ratio_phi+0.5*phi2-0.5*phi2star);
  double randnum = runif(1)[0];
  if(randnum<ratio_phi){
    at_phi += 1;
    for(int i=0; i<ND; ++i){
      phi_out[i]=phi_star(i,0);
    }
  }else{
    phi_out=miu_phi;
  }
  return phi_out;
}


// [[Rcpp::export]]
NumericVector f_zeta(NumericVector tran, NumericVector zeta, NumericVector s, NumericVector at_zeta, double c_zeta, int N, int NT, int NS){
  NumericVector zeta_out(NS*NS);
  NumericVector s1(N);
  NumericVector pri_zeta(NS);
  pri_zeta(0) = 1;
  pri_zeta(1) = 0;
  
  for(int i=0; i<NS; ++i){
    for(int j=0; j<(NS-1); ++j){
      double temp_zeta=0;
      int n_zeta =0;
      double miu_zeta=zeta[i*NS+j];
      for(int k=j; k<NS; ++k){
        for(int l=1; l<NT; ++l){
          for(int m=0; m<N; ++m){
            if((s[l*N+m]==k)&&(s[(l-1)*N+m]==i)){
              temp_zeta += exp(tran[l*N*NS*NS+m*NS*NS]-zeta[0])/pow((1+exp(tran[l*N*NS*NS+m*NS*NS]-zeta[0])),2);
              if(k==j){
                n_zeta +=1;
              }
            }
          }
        }
      }
      double sigma_zeta = 1/(temp_zeta+1);
      double zeta_star = rnorm(1,miu_zeta,sqrt(c_zeta*sigma_zeta))[0];
      
      double ratio1=0;
      double ratio2=0;
      for(int k=j; k<NS; ++k){
        for(int l=1; l<NT; ++l){
          for(int m=0; m<N; ++m){
            if((s[l*N+m]==k)&&(s[(l-1)*N+m]==i)){
              ratio1 += log(1+exp(tran[l*N*NS*NS+m*NS*NS]-zeta[0]+miu_zeta));
              ratio2 += log(1+exp(tran[l*N*NS*NS+m*NS*NS]-zeta[0]+zeta_star));
            }
          }
        }
      }
      
      double ratio_zeta = exp((zeta_star-miu_zeta)*n_zeta+ratio1-ratio2+0.5*pow(miu_zeta-pri_zeta(i),2)-0.5*pow(zeta_star-pri_zeta(i),2));
      double randnum = runif(1)[0];
      if(randnum<ratio_zeta){
        at_zeta[i*NS+j] += 1;
        zeta_out[i*NS+j]=zeta_star;
      }else{
        zeta_out[i*NS+j]=miu_zeta;
      }
    }
  }
  return zeta_out;
}


// [[Rcpp::export]]
NumericVector f_Bbeta(NumericVector fc, NumericVector alpha, int N, int NT, int NS, int NF){
  
  NumericVector Bbeta_out(NS*N*NT);
  for(int i=0; i<NS; ++i){
    for(int j=0; j<(N*NT); ++j){
      NumericVector fc_temp(NF+1);
      NumericVector alpha_temp(NF+1);
      for(int h=0; h<(NF+1); ++h){
        fc_temp[h]=fc[h*N*NT+j];
        // alpha_temp[h]=alpha[i*(NF+1)+h];
        alpha_temp[h]=alpha[h*NS+i];
      }
      Bbeta_out[i*N*NT+j]+=sum(fc_temp*alpha_temp);
    }
  }
  return(Bbeta_out);
}


// [[Rcpp::export]]
List f_alpha(NumericVector alpha, NumericVector fc, NumericVector y, NumericVector Bbeta, 
             NumericVector s, NumericVector sigma, NumericVector betan, NumericVector lambda2, 
             NumericVector tau2, NumericVector pi0, NumericVector pi1, NumericVector Z, 
             NumericVector sg, NumericVector groups, NumericVector betaprior, int N, int NT, 
             int NS, int NF, int ng, int eg){
  
  NumericVector alpha_out(NS*(eg * ng + 1));
  NumericVector betan_out(NS*ng);
  NumericVector tau2_out(NS*ng);
  NumericVector pi1_out(NS*ng);
  NumericVector Z_out(NS*ng);
  NumericVector lambda2_out(NS*ng);
  NumericVector pi0_out(NS);
  for(int i=0; i<NS; ++i){
    NumericVector alpha_used(eg * ng);
    double priorb = betaprior(i);
    double intercept_used;
    intercept_used = alpha((eg*ng)*NS + i);
    double pi0_used = pi0(i);
    double sig2_used = sigma(i);
    NumericVector Z_used(ng);
    NumericVector betan_used(ng);
    NumericVector tau2_used(ng);
    NumericVector pi1_used(ng);
    NumericVector lambda2_used(ng);
    int num = table(s)[i];
    NumericVector y_used(num);
    NumericVector fc_used(num*(eg*ng + 1));
    NumericMatrix sub_x(num, eg * ng);
    NumericMatrix temp_fc1(eg*ng+1,num); // used for intercept
    NumericMatrix temp_fc2(num,(eg*ng)+1);
    NumericMatrix temp(num,1);
    for(int q=0; q<(eg*ng); ++q){
      alpha_used(q) = alpha(q*NS + i); // every iter for ig, need to give value repeatedly
    }
    for(int p=0; p<ng; ++p){
      Z_used(p) = Z(p*NS + i);
      betan_used(p) = betan(p*NS + i); // every iter for ig, need to give value repeatedly
      tau2_used(p) = tau2(p*NS + i);
      pi1_used(p) = pi1(p*NS + i);
      lambda2_used(p) = lambda2(p*NS + i);
    }
    
    int cal = 0;
    for(int j=0; j<(N*NT); ++j) {
      if(s[j]==i){
        y_used(cal) = y(j);
        for(int k=0; k<(eg*ng+1); ++k){
          temp_fc1(k,cal) = fc[k*N*NT+j];
          temp_fc2(cal,k) = fc[k*N*NT+j];
        }
        temp(cal,0) = y[j]-Bbeta[i*N*NT+j];
        for(int gg=0; gg<(ng*eg + 1); ++gg){
          fc_used(gg*num + cal) = fc(gg*N*NT+j);
        }
        for(int l=0; l<(eg * ng); ++l){
          sub_x(cal,l) = fc(l*N*NT+j);
        }
        cal += 1;
      } 
    }
    
    
    for (int ig = 0; ig < ng; ++ig) {
      NumericMatrix sub_x1(num, eg); // represent x that we use
      NumericMatrix sub_x2(eg, num); // represent transpose of sub_x1
      NumericVector ind(eg); // represent col index
      for(int ee = 0; ee < eg; ++ee){
        ind(ee) = groups(ig*eg + ee) - 1; // -1 for index
      }
      
      for(int nt = 0; nt < num; ++nt){
        for(int k=0; k<eg; ++k){
          sub_x1(nt,k) = fc_used(ind(k)*num + nt);
          sub_x2(k,nt) = fc_used(ind(k)*num + nt);
        }
      }
      
      NumericMatrix xx(eg, eg);
      xx = f_mmult(sub_x2, sub_x1);
      double tau_used = tau2_used(ig);
      NumericMatrix D(eg, eg);
      NumericMatrix A(eg, eg);
      NumericMatrix AA(eg, eg);
      for(int m=0; m<eg; ++m){
        D(m, m) = 1 / tau_used;
      }
      A = f_madd(xx, D);
      AA = f_madd(A, f_transpose(A)) / 2;
      for(int ee = 0; ee < eg; ++ee){
        alpha_used(ig * eg + ee) = 0;
      }
      Z_used(ig) = 0;
      betan_used(ig) = 0;
      NumericMatrix b(num, 1);
      NumericMatrix xb(eg, 1);
      NumericVector subsub_x(eg * ng);
      for(int cc = 0; cc < num; ++cc){
        for(int l=0; l<(eg * ng); ++l){
          subsub_x(l) = sub_x(cc,l);
        }
        b(cc, 0) = y_used(cc) - intercept_used - sum(subsub_x * alpha_used);
      }
      xb = f_mmult(sub_x2, b);
      NumericMatrix AA_inv = f_inv(AA);
      NumericMatrix AAxb(eg, 1);
      AAxb = f_mmult(AA_inv, xb);
      NumericMatrix xbAAxb = f_mmult(f_transpose(xb), AAxb);
      double xbAAxbvalue = xbAAxb(0, 0);
      double maxAA = f_maxValue(AA);
      double L1 = (-sg(ig) / 2) * log(tau_used);
      double L2 = (-0.5) * log(f_det(AA));
      double L3 = (0.5 / sig2_used) * xbAAxbvalue;
      double L = (-sg(ig) / 2) * log(tau_used) + (-0.5) * log(f_det(AA / maxAA)) + (-eg / 2) * log(maxAA) + (0.5 / sig2_used) * xbAAxbvalue;
      double pi1_cal = pi0_used / (pi0_used + (1 - pi0_used) * exp(L));
      double randnum = runif(1)[0];
      NumericMatrix alpha_star(eg, 1);
      NumericVector alpha_star1(eg);
      
      double new_tau2;
      if (randnum >= pi1_cal){
        alpha_star = f_mvnorm(1, AAxb, sig2_used * AA_inv);
        for(int ee = 0; ee < eg; ++ee){
          alpha_star1(ee) = alpha_star(ee, 0);
        }
        Z_used(ig) = 1;
        betan_used(ig) = sum(alpha_star1 * alpha_star1);
        double a1 = sqrt(lambda2_used(ig) * sig2_used / sum(alpha_star1 * alpha_star1));
        double a2 = lambda2_used(ig);
        double new_tau2 = 1 / (f_invgauss(1, a1, a2)[0]) + 1e-10;
        tau2_used(ig) = new_tau2;
      } else {
        double shape = (sg(ig) + 1.0) / 2.0;
        double rate = lambda2_used(ig) / 2.0;
        double new_tau2 = 1 / (rgamma(1, shape, 1 / rate)[0]) + 1e-10;
        tau2_used(ig) = new_tau2;
      }
      for(int ee = 0; ee < eg; ++ee){
        alpha_used(ee + eg * ig) = alpha_star(ee, 0);
      }
      pi1_used(ig) = pi1_cal;
      
      double paral1 = 4 + ((sg(ig) + 1) / 2);
      double paral2 = 1 + (new_tau2 / 2);
      double new_lambda = rgamma(1, paral1, 1 / paral2)[0];
      lambda2_used(ig) = new_lambda;
    }
    for(int en = 0; en < (eg * ng); ++en){
      alpha_out(en*NS + i) = alpha_used(en);
    }
    
    // for intercept
    double mu_post = sum(y_used) / (sigma(i) + num);
    double sig_post = sigma(i) / (sigma(i) + num);
    double new_intercept = rnorm(1, mu_post, sqrt(sig_post))[0];
    alpha_out[(eg*ng)*NS + i]=new_intercept;
    
    
    for(int ig = 0; ig < ng; ++ig){
      betan_out(ig*NS + i) = betan_used(ig);
      tau2_out(ig*NS + i) = tau2_used(ig);
      pi1_out(ig*NS + i) = pi1_used(ig);
      Z_out(ig*NS + i) = Z_used(ig);
      lambda2_out(ig*NS + i) = lambda2_used(ig);
    }
    
    double new_pi0;
    double paraA = 10 + ng - sum(Z_used);
    double paraB = 1 + sum(Z_used);
    new_pi0 = rbeta(1, paraA, paraB)[0];
    pi0_out[i] = new_pi0;
  }
  
  return List::create(Rcpp::Named("new_alpha") = alpha_out,
                      Rcpp::Named("new_betan") = betan_out,
                      Rcpp::Named("new_tau") = tau2_out,
                      Rcpp::Named("new_pi1") = pi1_out,
                      Rcpp::Named("new_Z") = Z_out,
                      Rcpp::Named("new_lambda") = lambda2_out,
                      Rcpp::Named("new_pi0") = pi0_out); 
}


// [[Rcpp::export]]
NumericVector f_sigma(NumericVector y, NumericVector Bbeta, NumericVector betan, NumericVector tau2,
                      NumericVector Z, NumericVector sg, NumericVector s, int N, int NT, int NS, int ng) {
  NumericVector sigma_out(NS);
  NumericVector sigma_sum(NS);
  NumericVector Z_used(ng);
  NumericVector betan_used(ng);
  NumericVector tau2_used(ng);
  NumericVector betan_tau(ng);
  for(int i=0; i<NS; ++i){
    for(int p=0; p<ng; ++p){
      Z_used(p) = Z(p*NS + i);
      betan_used(p) = betan(p*NS + i); 
      tau2_used(p) = tau2(p*NS + i); 
    }
    for(int ig = 0; ig < ng; ++ig){
      betan_tau(ig) = betan_used(ig) / tau2_used(ig);
      // betan_tau(ig) = (betan_used(ig) * betan_used(ig))/ tau2_used(ig);
    }
    for(int j=0; j<(N*NT); ++j){
      if(s[j]==i){
        sigma_sum[i] += pow((y[j]-Bbeta[i*N*NT+j]),2);
      }
    }
    sigma_out[i]=1/(rgamma(1,9+0.5*table(s)[i]+0.5*sum(Z_used*sg),1/(4+0.5*sigma_sum[i]+0.5*sum(betan_tau)))[0]);
    // sigma_out[i]=1/(rgamma(1,9+0.5*table(s)[i],1/(4+0.5*sigma_sum[i]))[0]);
  }
  return sigma_out;  
}


// [[Rcpp::export]]
List f_label(NumericVector s, NumericVector sigma, NumericVector alpha, NumericVector betan, 
             NumericVector tau2, NumericVector pi1, NumericVector Z, NumericVector lambda2, 
             NumericVector pi0, int N, int NT, int NF, int NS, int eg, int ng){
  NumericVector s_label(N*NT);
  NumericVector sigma_label(NS);
  NumericVector alpha_label(NS*(eg * ng + 1));
  NumericVector betan_label(NS*ng);
  NumericVector tau2_label(NS*ng);
  NumericVector pi1_label(NS*ng);
  NumericVector Z_label(NS*ng);
  NumericVector lambda2_label(NS*ng);
  NumericVector pi0_label(NS);
  NumericVector test_miu(NS);
  NumericVector test_order(NS);
  NumericVector test_order1(NS);
  NumericVector test(NS);
  for(int i=0; i<NS; ++i){
    test_miu[i]=alpha[(eg*ng)*NS + i];
    test[i]=i+1;
  }
  for(int i=0; i<NS; ++i){
    test_order[i]=f_order(test_miu)(i,0);
  }
  for(int i=0; i<NS; ++i){
    test_order1[i]=f_order(test_order)(i,0);  
  }
  if(is_true(any(test_order!=test))){
    for(int i=0; i<(N*NT); ++i){
      s_label[i]=test_order1[s[i]]-1;
    }
    for(int i=0; i<NS; ++i){
      sigma_label[i]=sigma[test_order[i]-1];
      pi0_label[i]=pi0[test_order[i]-1];
    }
    for(int i=0; i<NS; ++i){
      for(int j=0; j<(ng*eg+1); ++j){
        alpha_label[j*NS+i]=alpha[j*NS+(test_order[i]-1)];
      }
      for(int k=0; k<ng; ++k){
        betan_label[k*NS+i]=betan[k*NS+(test_order[i]-1)];
        tau2_label[k*NS+i]=tau2[k*NS+(test_order[i]-1)];
        pi1_label[k*NS+i]=pi1[k*NS+(test_order[i]-1)];
        Z_label[k*NS+i]=Z[k*NS+(test_order[i]-1)];
        lambda2_label[k*NS+i]=lambda2[k*NS+(test_order[i]-1)];
      }
    }
  }else{
    s_label=s;
    sigma_label=sigma;
    alpha_label=alpha;
    betan_label=betan;
    tau2_label=tau2;
    pi1_label=pi1;
    Z_label=Z;
    lambda2_label=lambda2;
    pi0_label=pi0;
  }
  return List::create(Rcpp::Named("s_label") = s_label,
                      Rcpp::Named("sigma_label") = sigma_label,
                      Rcpp::Named("alpha_label") = alpha_label,
                      Rcpp::Named("betan_label") = betan_label,
                      Rcpp::Named("tau2_label") = tau2_label,
                      Rcpp::Named("pi1_label") = pi1_label,
                      Rcpp::Named("Z_label") = Z_label,
                      Rcpp::Named("lambda2_label") = lambda2_label,
                      Rcpp::Named("pi0_label") = pi0_label);
}



// [[Rcpp::export]]
List mcmc(NumericVector y, NumericVector x, NumericVector fc, NumericVector d, NumericVector alpha,
          NumericVector sigma, NumericVector tau, NumericVector zeta, NumericVector phi, 
          NumericVector betan, NumericVector lambda2, NumericVector tau2, NumericVector pi0, 
          NumericVector pi1, NumericVector Z, NumericVector sg, NumericVector groups, NumericVector s,
          NumericVector betaprior, int N, int NT, int NS, int NF, int ND, double c_tau, double c_zeta, 
          double c_phi, int iter, int rep, int T, int S, int ng, int eg){
  NumericVector Bbeta(NS*N*NT);
  NumericVector L(NS*N*NT);
  NumericVector tran(NT*N*NS*NS);
  NumericVector p(NT*N*NS*NS);
  NumericVector p0(NS);
  NumericVector at_tau(NS-1);
  NumericVector at_zeta(NS*NS);
  int at_phi = 0;
  
  NumericVector phi_str(S*ND);
  NumericVector zeta_str(S*NS*NS);
  NumericVector sigma_str(S*NS);
  NumericVector tau_str(S*(NS-1));
  NumericVector alpha_str(S*NS*(ng*eg+1));
  NumericVector s_str(S*N*NT);
  NumericVector pi1_str(S*NS*ng);
  
  for(int t=0; t<T; ++t){
    Bbeta = f_Bbeta(fc, alpha, N, NT, NS, NF);
    L = f_likeli(Bbeta, y, sigma, N, NT, NS);
    tran = f_tran(d, zeta, phi, N, NT, NS, ND);
    
    List prob = f_p(tran,tau,N,NT,NS);
    p = prob[0];
    p0 = prob[1];
    
    if(t>50){
      s = f_s(L, p, p0, N, NT, NS);
    }
    // Rcpp::Rcout << "The state vector is: " << s << std::endl;
    
    sigma = f_sigma(y, Bbeta, betan, tau2, Z, sg, s, N, NT, NS, ng);
    tau = f_tau(tau, s, at_tau, c_tau, N, NS);
    zeta = f_zeta(tran,zeta,s,at_zeta,c_zeta,N,NT,NS);
    phi = f_phi(d,zeta,phi,s,at_phi,c_phi,N,NT,ND,NS);
    Bbeta = f_Bbeta(fc, alpha, N, NT, NS, NF);
    List reg_result = f_alpha(alpha, fc, y, Bbeta, s, sigma, betan, lambda2, tau2, pi0, pi1, Z, sg, groups, betaprior, N, NT, NS, NF, ng, eg);
    alpha = reg_result[0];
    betan = reg_result[1];
    tau2 = reg_result[2];
    pi1 = reg_result[3];
    Z = reg_result[4];
    lambda2 = reg_result[5];
    pi0 = reg_result[6];
    
    List label=f_label(s, sigma, alpha, betan, tau2, pi1, Z, lambda2, pi0, N, NT, NF, NS, eg, ng);
    s=label[0];
    sigma=label[1];
    alpha=label[2];
    betan=label[3];
    tau2=label[4];
    pi1=label[5];
    Z=label[6];
    lambda2=label[7];
    pi0=label[8];
    
    if(t>(T-S-1)){
      int ind = t-(T-S);
      for(int i=0; i<(NS-1); ++i){
        tau_str[ind*(NS-1)+i] = tau[i];
      }
      for(int i=0; i<NS; ++i){
        for(int j=0; j<NS; ++j){
          zeta_str[ind*NS*NS+i*NS+j] = zeta[i*NS+j];
        }
      }
      for(int i=0; i<ND; ++i){
        phi_str[ind*ND+i] = phi[i];
      }
      for(int i=0; i<NT; ++i){
        for(int j=0; j<N; ++j){
          s_str[ind*N*NT+i*N+j] = s[i*N+j];
        }
      }
      for(int i=0; i<NS; ++i){
        for(int j=0; j<(NF+1); ++j){
          alpha_str[ind*NS*(NF+1)+i*(NF+1)+j] = alpha[j*NS + i];
        }
      }
      for(int i=0; i<NS; ++i){
        for(int j=0; j<ng; ++j){
          pi1_str[ind*NS*ng+i*ng+j] = pi1[j*NS + i];
        }
      }
      for(int i=0; i<NS; ++i){
        sigma_str[ind*NS+i] = sigma[i];
      }
    }
    
    
    if((t+1)%10==0){
      Rprintf("\rRunning MCMC, %d/%d replication, %f%% completed...",rep,iter,(t+1)*100.0/T);
    }
  }
  
  return List::create(Rcpp::Named("tau_str") = tau_str,
                      Rcpp::Named("zeta_str") = zeta_str,
                      Rcpp::Named("phi_str") = phi_str,
                      Rcpp::Named("s_str") = s_str,
                      Rcpp::Named("alpha_str") = alpha_str,
                      Rcpp::Named("pi1_str") = pi1_str,
                      Rcpp::Named("sigma_str") = sigma_str,
                      Rcpp::Named("at_tau") = at_tau,
                      Rcpp::Named("at_zeta") = at_zeta,
                      Rcpp::Named("at_phi") = at_phi);
}


// [[Rcpp::export]]
List mcmc_hmm(NumericVector y, NumericVector x, NumericVector fc, NumericVector alpha,
          NumericVector sigma, NumericVector betan, NumericVector lambda2, NumericVector tau2, 
          NumericVector pi0, NumericVector pi1, NumericVector Z, NumericVector sg, NumericVector groups, 
          NumericVector s, NumericVector betaprior, int N, int NT, int NS, int NF, int iter, int rep, 
          int T, int S, int ng, int eg, NumericVector ele){
  NumericVector Bbeta(NS*N*NT);
  NumericVector L(NS*N*NT);
  NumericVector p(NS*NS);
  NumericVector p0(NS);
  
  NumericVector p_str(S*NS*NS);
  NumericVector sigma_str(S*NS);
  NumericVector alpha_str(S*NS*(ng*eg+1));
  NumericVector s_str(S*N*NT);
  NumericVector pi1_str(S*NS*ng);
  
  NumericVector ones(NS);
  for(int i=0; i<NS; ++i){
    ones(i) = 1;
  }
  for(int t=0; t<T; ++t){
    Bbeta = f_Bbeta(fc, alpha, N, NT, NS, NF);
    L = f_likeli(Bbeta, y, sigma, N, NT, NS);
    p0 = rdirichlet(1, ones);
    // Rcpp::Rcout << "The t vector is: " << t << std::endl;
    // Rcpp::Rcout << "The p0 vector is: " << p0 << std::endl;
    p = f_trprob(s,N,NT,NS,ele);
    // Rcpp::Rcout << "The p vector is: " << p << std::endl;
    if(t>50){
      s = f_s_hmm(L, p, p0, N, NT, NS);
    }
    // Rcpp::Rcout << "The state vector is: " << s << std::endl;
    
    sigma = f_sigma(y, Bbeta, betan, tau2, Z, sg, s, N, NT, NS, ng);
    Bbeta = f_Bbeta(fc, alpha, N, NT, NS, NF);
    List reg_result = f_alpha(alpha, fc, y, Bbeta, s, sigma, betan, lambda2, tau2, pi0, pi1, Z, sg, groups, betaprior, N, NT, NS, NF, ng, eg);
    alpha = reg_result[0];
    betan = reg_result[1];
    tau2 = reg_result[2];
    pi1 = reg_result[3];
    Z = reg_result[4];
    lambda2 = reg_result[5];
    pi0 = reg_result[6];
    
    List label=f_label(s, sigma, alpha, betan, tau2, pi1, Z, lambda2, pi0, N, NT, NF, NS, eg, ng);
    s=label[0];
    sigma=label[1];
    alpha=label[2];
    betan=label[3];
    tau2=label[4];
    pi1=label[5];
    Z=label[6];
    lambda2=label[7];
    pi0=label[8];
    
    if(t>(T-S-1)){
      int ind = t-(T-S);
      for(int i=0; i<NS; ++i){
        for(int j=0; j<NS; ++j){
          p_str[ind*NS*NS+i*NS+j] = p[j*NS+i];
        }
      }
      for(int i=0; i<NT; ++i){
        for(int j=0; j<N; ++j){
          s_str[ind*N*NT+i*N+j] = s[i*N+j];
        }
      }
      for(int i=0; i<NS; ++i){
        for(int j=0; j<(NF+1); ++j){
          alpha_str[ind*NS*(NF+1)+i*(NF+1)+j] = alpha[j*NS + i];
        }
      }
      for(int i=0; i<NS; ++i){
        for(int j=0; j<ng; ++j){
          pi1_str[ind*NS*ng+i*ng+j] = pi1[j*NS + i];
        }
      }
      for(int i=0; i<NS; ++i){
        sigma_str[ind*NS+i] = sigma[i];
      }
    }
    
    
    if((t+1)%10==0){
      Rprintf("\rRunning MCMC, %d/%d replication, %f%% completed...",rep,iter,(t+1)*100.0/T);
    }
  }
  
  return List::create(Rcpp::Named("p_str") = p_str,
                      Rcpp::Named("s_str") = s_str,
                      Rcpp::Named("alpha_str") = alpha_str,
                      Rcpp::Named("pi1_str") = pi1_str,
                      Rcpp::Named("sigma_str") = sigma_str);
}

