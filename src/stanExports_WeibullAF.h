// Generated by rstantools.  Do not edit by hand.

/*
    survHE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    survHE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with survHE.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_WeibullAF_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_WeibullAF");
    reader.add_event(67, 65, "end", "model_WeibullAF");
    return reader;
}
template <typename T0__, typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
log_h(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
          const T1__& shape,
          const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 6;
        validate_non_negative_index("log_h", "num_elements(t)", num_elements(t));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> log_h(num_elements(t));
        stan::math::initialize(log_h, DUMMY_VAR__);
        stan::math::fill(log_h, DUMMY_VAR__);
        current_statement_begin__ = 7;
        stan::math::assign(log_h, subtract(add(stan::math::log(shape), multiply((shape - 1), stan::math::log(elt_divide(t, scale)))), stan::math::log(scale)));
        current_statement_begin__ = 8;
        return stan::math::promote_scalar<fun_return_scalar_t__>(log_h);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct log_h_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
          const T1__& shape,
          const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) const {
        return log_h(t, shape, scale, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
log_S(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
          const T1__& shape,
          const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 13;
        validate_non_negative_index("log_S", "num_elements(t)", num_elements(t));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> log_S(num_elements(t));
        stan::math::initialize(log_S, DUMMY_VAR__);
        stan::math::fill(log_S, DUMMY_VAR__);
        current_statement_begin__ = 14;
        for (int i = 1; i <= num_elements(t); ++i) {
            current_statement_begin__ = 15;
            stan::model::assign(log_S, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        -(pow((get_base1(t, i, "t", 1) / get_base1(scale, i, "scale", 1)), shape)), 
                        "assigning variable log_S");
        }
        current_statement_begin__ = 17;
        return stan::math::promote_scalar<fun_return_scalar_t__>(log_S);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct log_S_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
          const T1__& shape,
          const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) const {
        return log_S(t, shape, scale, pstream__);
    }
};
template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
surv_weibullAF_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
                        const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& d,
                        const T2__& shape,
                        const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 22;
        validate_non_negative_index("log_lik", "num_elements(t)", num_elements(t));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> log_lik(num_elements(t));
        stan::math::initialize(log_lik, DUMMY_VAR__);
        stan::math::fill(log_lik, DUMMY_VAR__);
        current_statement_begin__ = 23;
        local_scalar_t__ prob(DUMMY_VAR__);
        (void) prob;  // dummy to suppress unused var warning
        stan::math::initialize(prob, DUMMY_VAR__);
        stan::math::fill(prob, DUMMY_VAR__);
        current_statement_begin__ = 24;
        stan::math::assign(log_lik, add(elt_multiply(d, log_h(t, shape, scale, pstream__)), log_S(t, shape, scale, pstream__)));
        current_statement_begin__ = 25;
        stan::math::assign(prob, sum(log_lik));
        current_statement_begin__ = 26;
        return stan::math::promote_scalar<fun_return_scalar_t__>(prob);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
surv_weibullAF_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
                        const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& d,
                        const T2__& shape,
                        const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) {
    return surv_weibullAF_lpdf<false>(t,d,shape,scale, pstream__);
}
struct surv_weibullAF_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
                        const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& d,
                        const T2__& shape,
                        const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& scale, std::ostream* pstream__) const {
        return surv_weibullAF_lpdf(t, d, shape, scale, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_WeibullAF
  : public stan::model::model_base_crtp<model_WeibullAF> {
private:
        int n;
        vector_d t;
        vector_d d;
        int H;
        matrix_d X;
        vector_d mu_beta;
        vector_d sigma_beta;
        double a_alpha;
        double b_alpha;
public:
    model_WeibullAF(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_WeibullAF(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_WeibullAF_namespace::model_WeibullAF";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 31;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            current_statement_begin__ = 32;
            validate_non_negative_index("t", "n", n);
            context__.validate_dims("data initialization", "t", "vector_d", context__.to_vec(n));
            t = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("t");
            pos__ = 0;
            size_t t_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < t_j_1_max__; ++j_1__) {
                t(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 33;
            validate_non_negative_index("d", "n", n);
            context__.validate_dims("data initialization", "d", "vector_d", context__.to_vec(n));
            d = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("d");
            pos__ = 0;
            size_t d_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < d_j_1_max__; ++j_1__) {
                d(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 34;
            context__.validate_dims("data initialization", "H", "int", context__.to_vec());
            H = int(0);
            vals_i__ = context__.vals_i("H");
            pos__ = 0;
            H = vals_i__[pos__++];
            current_statement_begin__ = 35;
            validate_non_negative_index("X", "n", n);
            validate_non_negative_index("X", "H", H);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(n,H));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, H);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = H;
            size_t X_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 36;
            validate_non_negative_index("mu_beta", "H", H);
            context__.validate_dims("data initialization", "mu_beta", "vector_d", context__.to_vec(H));
            mu_beta = Eigen::Matrix<double, Eigen::Dynamic, 1>(H);
            vals_r__ = context__.vals_r("mu_beta");
            pos__ = 0;
            size_t mu_beta_j_1_max__ = H;
            for (size_t j_1__ = 0; j_1__ < mu_beta_j_1_max__; ++j_1__) {
                mu_beta(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 37;
            validate_non_negative_index("sigma_beta", "H", H);
            context__.validate_dims("data initialization", "sigma_beta", "vector_d", context__.to_vec(H));
            sigma_beta = Eigen::Matrix<double, Eigen::Dynamic, 1>(H);
            vals_r__ = context__.vals_r("sigma_beta");
            pos__ = 0;
            size_t sigma_beta_j_1_max__ = H;
            for (size_t j_1__ = 0; j_1__ < sigma_beta_j_1_max__; ++j_1__) {
                sigma_beta(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "sigma_beta", sigma_beta, 0);
            current_statement_begin__ = 38;
            context__.validate_dims("data initialization", "a_alpha", "double", context__.to_vec());
            a_alpha = double(0);
            vals_r__ = context__.vals_r("a_alpha");
            pos__ = 0;
            a_alpha = vals_r__[pos__++];
            check_greater_or_equal(function__, "a_alpha", a_alpha, 0);
            current_statement_begin__ = 39;
            context__.validate_dims("data initialization", "b_alpha", "double", context__.to_vec());
            b_alpha = double(0);
            vals_r__ = context__.vals_r("b_alpha");
            pos__ = 0;
            b_alpha = vals_r__[pos__++];
            check_greater_or_equal(function__, "b_alpha", b_alpha, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 43;
            validate_non_negative_index("beta", "H", H);
            num_params_r__ += H;
            current_statement_begin__ = 44;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_WeibullAF() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 43;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "H", H);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(H));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(H);
        size_t beta_j_1_max__ = H;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 44;
        if (!(context__.contains_r("alpha")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha", "double", context__.to_vec());
        double alpha(0);
        alpha = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, alpha);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 43;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(H, lp__);
            else
                beta = in__.vector_constrain(H);
            current_statement_begin__ = 44;
            local_scalar_t__ alpha;
            (void) alpha;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha = in__.scalar_lb_constrain(0, lp__);
            else
                alpha = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 48;
            validate_non_negative_index("linpred", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> linpred(n);
            stan::math::initialize(linpred, DUMMY_VAR__);
            stan::math::fill(linpred, DUMMY_VAR__);
            current_statement_begin__ = 49;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 50;
            stan::math::assign(linpred, multiply(X, beta));
            current_statement_begin__ = 51;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 52;
                stan::model::assign(mu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::exp(get_base1(linpred, i, "linpred", 1)), 
                            "assigning variable mu");
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 48;
            size_t linpred_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < linpred_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(linpred(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: linpred" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable linpred: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 49;
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(mu(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: mu" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable mu: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 57;
            lp_accum__.add(gamma_log<propto__>(alpha, a_alpha, b_alpha));
            current_statement_begin__ = 58;
            lp_accum__.add(normal_log<propto__>(beta, mu_beta, sigma_beta));
            current_statement_begin__ = 59;
            lp_accum__.add(surv_weibullAF_lpdf<propto__>(t, d, alpha, mu, pstream__));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta");
        names__.push_back("alpha");
        names__.push_back("linpred");
        names__.push_back("mu");
        names__.push_back("scale");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(H);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_WeibullAF_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(H);
        size_t beta_j_1_max__ = H;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        double alpha = in__.scalar_lb_constrain(0);
        vars__.push_back(alpha);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 48;
            validate_non_negative_index("linpred", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> linpred(n);
            stan::math::initialize(linpred, DUMMY_VAR__);
            stan::math::fill(linpred, DUMMY_VAR__);
            current_statement_begin__ = 49;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<double, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 50;
            stan::math::assign(linpred, multiply(X, beta));
            current_statement_begin__ = 51;
            for (int i = 1; i <= n; ++i) {
                current_statement_begin__ = 52;
                stan::model::assign(mu, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::exp(get_base1(linpred, i, "linpred", 1)), 
                            "assigning variable mu");
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t linpred_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < linpred_j_1_max__; ++j_1__) {
                    vars__.push_back(linpred(j_1__));
                }
                size_t mu_j_1_max__ = n;
                for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                    vars__.push_back(mu(j_1__));
                }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 63;
            double scale;
            (void) scale;  // dummy to suppress unused var warning
            stan::math::initialize(scale, DUMMY_VAR__);
            stan::math::fill(scale, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 64;
            stan::math::assign(scale, stan::math::exp(get_base1(beta, 1, "beta", 1)));
            // validate, write generated quantities
            current_statement_begin__ = 63;
            vars__.push_back(scale);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_WeibullAF";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = H;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t linpred_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < linpred_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "linpred" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "scale";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = H;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t linpred_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < linpred_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "linpred" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t mu_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "mu" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "scale";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_WeibullAF_namespace::model_WeibullAF stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
