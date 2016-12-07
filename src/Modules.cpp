#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4Exponential_mod) {


    class_<rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> >("model_Exponential")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_Exponential_namespace::model_Exponential, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4Gamma_mod) {


    class_<rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> >("model_Gamma")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_Gamma_namespace::model_Gamma, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4GenF_mod) {


    class_<rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> >("model_GenF")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_GenF_namespace::model_GenF, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4GenGamma_mod) {


    class_<rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> >("model_GenGamma")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_GenGamma_namespace::model_GenGamma, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4Gompertz_mod) {


    class_<rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> >("model_Gompertz")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_Gompertz_namespace::model_Gompertz, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4PolyWeibull_mod) {


    class_<rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> >("model_PolyWeibull")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_PolyWeibull_namespace::model_PolyWeibull, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4RP_mod) {


    class_<rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> >("model_RP")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_RP_namespace::model_RP, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4WeibullAF_mod) {


    class_<rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> >("model_WeibullAF")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_WeibullAF_namespace::model_WeibullAF, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4WeibullPH_mod) {


    class_<rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> >("model_WeibullPH")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_WeibullPH_namespace::model_WeibullPH, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logLogistic_mod) {


    class_<rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> >("model_logLogistic")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logLogistic_namespace::model_logLogistic, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4logNormal_mod) {


    class_<rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> >("model_logNormal")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_logNormal_namespace::model_logNormal, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
