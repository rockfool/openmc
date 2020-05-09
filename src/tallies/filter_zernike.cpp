#include "openmc/tallies/filter_zernike.h"

#include <cmath>
#include <sstream>
#include <utility>  // For pair

#include <fmt/core.h>
#include <gsl/gsl>

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// ZernikeFilter implementation
//==============================================================================

void
ZernikeFilter::from_xml(pugi::xml_node node)
{
  set_order(std::stoi(get_node_value(node, "order")));
  x_ = std::stod(get_node_value(node, "x"));
  y_ = std::stod(get_node_value(node, "y"));
  r_ = std::stod(get_node_value(node, "r"));
}

void
ZernikeFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                            FilterMatch& match) const
{
  // Determine the normalized (r,theta) coordinates.
  double x = p->r().x - x_;
  double y = p->r().y - y_;
  double r = std::sqrt(x*x + y*y) / r_;
  double theta = std::atan2(y, x);

  if (r <= 1.0) {
    // Compute and return the Zernike weights.
    std::vector<double> zn(n_bins_);
    calc_zn(order_, r, theta, zn.data());
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i);
      match.weights_.push_back(zn[i]);
    }
  }
}

void
ZernikeFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  write_dataset(filter_group, "x", x_);
  write_dataset(filter_group, "y", y_);
  write_dataset(filter_group, "r", r_);
}

std::string
ZernikeFilter::text_label(int bin) const
{
  Expects(bin >= 0 && bin < n_bins_);
  for (int n = 0; n < order_+1; n++) {
    int last = (n + 1) * (n + 2) / 2;
    if (bin < last) {
      int first = last - (n + 1);
      int m = -n + (bin - first) * 2;
      return fmt::format("Zernike expansion, Z{},{}", n, m);
    }
  }
  UNREACHABLE();
}

void
ZernikeFilter::set_order(int order)
{
  if (order < 0) {
    throw std::invalid_argument{"Zernike order must be non-negative."};
  }
  order_ = order;
  n_bins_ = ((order+1) * (order+2)) / 2;
}

//==============================================================================
// ZernikeRadialFilter implementation
//==============================================================================

void
ZernikeRadialFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                                  FilterMatch& match) const
{
  // Determine the normalized radius coordinate.
  double x = p->r().x - x_;
  double y = p->r().y - y_;
  double r = std::sqrt(x*x + y*y) / r_;

  if (r <= 1.0) {
    // Compute and return the Zernike weights.
    std::vector<double> zn(n_bins_);
    calc_zn_rad(order_, r, zn.data());
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i);
      match.weights_.push_back(zn[i]);
    }
  }
}

std::string
ZernikeRadialFilter::text_label(int bin) const
{
  return "Zernike expansion, Z" + std::to_string(2*bin) + ",0";
}

void
ZernikeRadialFilter::set_order(int order)
{
  ZernikeFilter::set_order(order);
  n_bins_ = order / 2 + 1;
}


//==============================================================================
// MultipleZernikeFilter implementation
//==============================================================================

void
MultipleZernikeFilter::from_xml(pugi::xml_node node)
{
  auto orders = get_node_array<int>(node, "order");
  auto xs = get_node_array<double>(node, "x");
  auto ys = get_node_array<double>(node, "y");
  auto rs = get_node_array<double>(node, "r");
  len_ = orders.size();
  n_counter_ = -1;
  if(len_ == xs.size() && len_ == ys.size() && len_ == xs.size()){
    for (int i = 0; i < len_; i++){
      //orders_.push_back(orders[i]); 
      set_orders(orders);
      xs_.push_back(xs[i]);
      ys_.push_back(ys[i]);
      rs_.push_back(rs[i]);
    }
    update_n_bins(0);
  } else {
    throw std::invalid_argument{"MultipleZernike xs ys rs and orders \
                                 must be in same dimension."};
  }
}

void
MultipleZernikeFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                            FilterMatch& match, int n)
{
  // Determine the normalized (r,theta) coordinates.
  double xs = p->r().x - xs_[n];
  double ys = p->r().y - ys_[n];
  double rs = std::sqrt(xs*xs + ys*ys) / rs_[n];
  double theta = std::atan2(ys, xs);

  if (rs <= 1.0) {
    // Compute and return the Zernike weights.
    update_n_bins(n); 
    // update n_bins_
    std::vector<double> zn(n_bins_);
    calc_zn(orders_[n], rs, theta, zn.data());
    match.bins_.clear();
    match.weights_.clear();
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i);
      match.weights_.push_back(zn[i]);
    }
  }
}

void
MultipleZernikeFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                            FilterMatch& match) const  
{
  int n = n_counter_;
  // Determine the normalized (r,theta) coordinates.
  double xs = p->r().x - xs_[n];
  double ys = p->r().y - ys_[n];
  double rs = std::sqrt(xs*xs + ys*ys) / rs_[n];
  double theta = std::atan2(ys, xs);
  
  if (rs <= 1.0) {
    // Compute and return the Zernike weights.
    std::vector<double> zn(n_bins_);
    calc_zn(orders_[n], rs, theta, zn.data());
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i);
      match.weights_.push_back(zn[i]);
    }
  }
}

void
MultipleZernikeFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", orders_);
  write_dataset(filter_group, "x", xs_);
  write_dataset(filter_group, "y", ys_);
  write_dataset(filter_group, "r", rs_);
  write_dataset(filter_group, "len", len_);
}

std::string
MultipleZernikeFilter::text_label(int bin, int n) 
{
  update_n_bins(n);
  Expects(bin >= 0 && bin < n_bins_);
  for (int i = 0; i < orders_[n] + 1; i++) {
    int last = (i + 1) * (i + 2) / 2;
    if (bin < last) {
      int first = last - (i + 1);
      int m = -i + (bin - first) * 2;
      return fmt::format("Zernike expansion, Z{},{}", i, m);
    }
  }
  UNREACHABLE();
}

std::string
MultipleZernikeFilter::text_label(int bin) const
{
  int n = n_counter_;
  Expects(bin >= 0 && bin < n_bins_);
  for (int i = 0; i < orders_[n] + 1; i++) {
    int last = (i + 1) * (i + 2) / 2;
    if (bin < last) {
      int first = last - (i + 1);
      int m = -i + (bin - first) * 2;
      return fmt::format("Zernike expansion, Z{},{}", i, m);
    }
  }
  UNREACHABLE();
}

void
MultipleZernikeFilter::set_orders(gsl::span<const int32_t> orders)
{ 
  orders_.clear();
  v_n_bins_.clear();
  len_ = orders.size();
  n_counter_ = -1;
  for (int i = 0; i < len_; i++){
    if (orders[i] < 0) {
      throw std::invalid_argument{"MultipleZernike order must be non-negative."};
    }
    orders_.push_back(orders[i]);
    v_n_bins_.push_back(((orders[i] + 1) * (orders[i] + 2)) / 2);
  }
  update_n_bins(0);
}

void
MultipleZernikeFilter::set_xs(gsl::span<const double> xs)
{
  if (xs.size() == len_){
    xs_.clear();
    for (int i = 0; i < xs.size(); i++){
      xs_.push_back(xs[i]);
    }
  } else{
    throw std::invalid_argument{"MultipleZernike xs must match orders."};
  }
}

void
MultipleZernikeFilter::set_ys(gsl::span<const double> ys)
{
  if (ys.size() == len_){
    ys_.clear();
    for (int i = 0; i < ys.size(); i++){
      ys_.push_back(ys[i]);
    }
  } else {
    throw std::invalid_argument{"MultipleZernike ys must match orders.."};
  }
}

void
MultipleZernikeFilter::set_rs(gsl::span<const double> rs)
{
  if (rs.size() == len_){
    rs_.clear();
    for (int i = 0; i < rs.size(); i++){
      rs_.push_back(rs[i]);
    }
  } else {
    throw std::invalid_argument{"MultipleZernike rs must match orders."};
  }
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, ZernikeFilter*>
check_zernike_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<ZernikeFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a Zernike filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int
openmc_zernike_filter_get_order(int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int
openmc_zernike_filter_get_params(int32_t index, double* x, double* y,
                                 double* r)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the params.
  *x = filt->x();
  *y = filt->y();
  *r = filt->r();
  return 0;
}

extern "C" int
openmc_zernike_filter_set_order(int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->set_order(order);
  return 0;
}

extern "C" int
openmc_zernike_filter_set_params(int32_t index, const double* x,
                                 const double* y, const double* r)
{
  // Check the filter.
  auto check_result = check_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  if (x) filt->set_x(*x);
  if (y) filt->set_y(*y);
  if (r) filt->set_r(*r);
  return 0;
}

std::pair<int, MultipleZernikeFilter*>
check_multiple_zernike_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MultipleZernikeFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a Zernike filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int
openmc_multiple_zernike_filter_get_orders(int32_t index, int** orders, size_t* n)
{
  // Check the filter.
  auto check_result = check_multiple_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the order.
  *orders = filt->orders().data();
  *n = filt->len();
  return 0;
}

extern "C" int
openmc_multiple_zernike_filter_get_params(int32_t index, double** xs, double** ys,
                                 double** rs, size_t* n)
{
  // Check the filter.
  auto check_result = check_multiple_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the params.
  *xs = filt->xs().data();
  *ys = filt->ys().data();
  *rs = filt->rs().data();
  *n = filt->len();
  return 0;
}

extern "C" int
openmc_multiple_zernike_filter_set_orders(int32_t index, size_t n, const int* orders)
{
  // Check the filter.
  auto check_result = check_multiple_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->set_orders({orders, n});
  return 0;
}

extern "C" int
openmc_multiple_zernike_filter_set_params(int32_t index, size_t n, const double* xs,
                                 const double* ys, const double* rs)
{
  // Check the filter.
  auto check_result = check_multiple_zernike_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  if (xs) filt->set_xs({xs, n});
  if (ys) filt->set_ys({ys, n});
  if (rs) filt->set_rs({rs, n});
  return 0;
}

} // namespace openmc
