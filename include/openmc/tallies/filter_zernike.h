#ifndef OPENMC_TALLIES_FILTER_ZERNIKE_H
#define OPENMC_TALLIES_FILTER_ZERNIKE_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Gives Zernike polynomial moments of a particle's position
//==============================================================================

class ZernikeFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~ZernikeFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "zernike";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }
  virtual void set_order(int order);

  double x() const { return x_; }
  void set_x(double x) { x_ = x; }

  double y() const { return y_; }
  void set_y(double y) { y_ = y; }

  double r() const { return r_; }
  void set_r(double r) { r_ = r; }

  //----------------------------------------------------------------------------
  // Data members

protected:
  //! Cartesian x coordinate for the origin of this expansion.
  double x_;

  //! Cartesian y coordinate for the origin of this expansion.
  double y_;

  //! Maximum radius from the origin covered by this expansion.
  double r_;

  int order_;
};

class MultipleZernikeFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~MultipleZernikeFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "multiplezernike";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match, \
                    int n);
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin, int n) const override;

  //----------------------------------------------------------------------------
  // Accessors
  
  std::vector<int>& orders() const { return orders_; }
  virtual void set_orders(gsl::span<const int32_t> orders);

  std::vector<double>& xs() const { return xs_; }
  virtual void set_xs(gsl::span<const int32_t> xs);

  std::vector<double>& ys() const { return ys_; }
  virtual void set_ys(gsl::span<const int32_t> ys);

  std::vector<double>& rs() const { return rs_; }
  virtual void set_rs(gsl::span<const int32_t> rs);
  
  std::vector<int>& n_bins() const { return n_bins_; }
  
  int len() { return len_; }

  //----------------------------------------------------------------------------
  // Data members

protected:
  //! Cartesian x coordinate for the origin of this expansion.
  std::vector<double> xs_;

  //! Cartesian y coordinate for the origin of this expansion.
  std::vector<double> ys_;

  //! Maximum radius from the origin covered by this expansion.
  std::vector<double> rs_;
  
  //! Numbers of element in Cartesian coordinate
  int len_;
  
  //! Orders of expansion
  std::vector<int> orders_;
  
  //! override the n_bins_
  std::vector<int> n_bins_; 
};

//==============================================================================
//! Gives even order radial Zernike polynomial moments of a particle's position
//==============================================================================

class ZernikeRadialFilter : public ZernikeFilter
{
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "zernikeradial";}

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_order(int order) override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ZERNIKE_H
