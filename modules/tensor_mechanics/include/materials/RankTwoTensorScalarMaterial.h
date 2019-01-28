#ifndef RANKTWOTENSORSCALARMATERIAL_H
#define RANKTWOTENSORSCALARMATERIAL_H

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

class RankTwoTensorScalarMaterial;

template <> InputParameters validParams<RankTwoTensorScalarMaterial>();

class RankTwoTensorScalarMaterial : public Material {
public:
  RankTwoTensorScalarMaterial(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const std::string _new_material_property_name;

  const MaterialPropertyName _rank_two_tensor_name;
  const MaterialProperty<RankTwoTensor> &_input_rank_two_tensor;

  const bool _stateful_flag;

  const bool _compute_rate_flag;

  const MooseEnum _scalar_type;

  /// whether or not selected_qp has been set
  const bool _has_selected_qp;

  /// The std::vector will be evaluated at this quadpoint only if defined
  const unsigned int _selected_qp;

  const Point _point1;
  const Point _point2;
  Point _input_direction;

  const unsigned int _n_new_material_properties;
  const unsigned int _n_new_material_properties_old;
  std::vector<MaterialProperty<Real> *> _new_material_properties;
  std::vector<const MaterialProperty<Real> *> _new_material_properties_old;

private:
  unsigned int nNewMP(const bool old = false);
  unsigned int check_qp(const unsigned int /*_qp*/) const;
};

#endif // RANKTWOTENSORSCALARMATERIAL_H
