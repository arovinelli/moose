#include "MooseError.h"
#include "RankTwoScalarTools.h"
#include "RankTwoTensorScalarMaterial.h"

registerMooseObject("MooseApp", RankTwoTensorScalarMaterial);

template <> InputParameters validParams<RankTwoTensorScalarMaterial>() {
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>(
      "rank_two_tensor_name", "The rank two material tensor names");

  params.addRequiredParam<std::string>("new_material_property_name",
                                       "The rank two material tensor names");

  params.addParam<bool>("stateful_flag",
                        "a boolean vector to decleare the new material "
                        "property as stateful, "
                        "default is false");
  params.addParam<bool>(
      "compute_rate_flag",
      "if we want to compute the rate of the selected scalar, default is "
      "false");
  params.addParam<MooseEnum>("scalar_type", RankTwoScalarTools::scalarOptions(),
                             "Type of scalar output");
  params.addParam<unsigned int>(
      "selected_qp", "Evaluate the tensor at this quadpoint.  This option only "
                     "needs to be used if "
                     "you are interested in a particular quadpoint in each "
                     "element: otherwise do "
                     "not include this parameter in your input file");
  params.addParamNamesToGroup("selected_qp", "Advanced");

  params.addParam<Point>("point1", Point(0, 0, 0),
                         "Start point for axis used to calculate some "
                         "cylinderical material tensor quantities");
  params.addParam<Point>(
      "point2", Point(0, 1, 0),
      "End point for axis used to calculate some material tensor quantities");
  params.addParam<Point>("direction", Point(0, 0, 1), "Direction vector");
  params.addClassDescription(
      "this material class add scalar material compute from a rank two tensor, "
      "like the stress using RankTwoScalarTools");

  return params;
}

RankTwoTensorScalarMaterial::RankTwoTensorScalarMaterial(
    const InputParameters &parameters)
    : Material(parameters),

      _new_material_property_name(
          parameters.get<std::string>("new_material_property_name")),

      _rank_two_tensor_name(
          parameters.get<MaterialPropertyName>("rank_two_tensor_name")),
      _input_rank_two_tensor(
          getMaterialPropertyByName<RankTwoTensor>(_rank_two_tensor_name)),
      _stateful_flag(parameters.isParamSetByUser("stateful_flag")
                         ? parameters.get<bool>("stateful_flag")
                         : false),

      _compute_rate_flag(parameters.isParamSetByUser("compute_rate_flag")
                             ? parameters.get<bool>("compute_rate_flag")
                             : false),

      _scalar_type(parameters.get<MooseEnum>("scalar_type")),
      _has_selected_qp(parameters.isParamSetByUser("selected_qp")),
      _selected_qp(_has_selected_qp ? getParam<unsigned int>("selected_qp")
                                    : 0),
      _point1(parameters.get<Point>("point1")),
      _point2(parameters.get<Point>("point2")),
      _input_direction(parameters.get<Point>("direction") /
                       parameters.get<Point>("direction").norm()),
      _n_new_material_properties(RankTwoTensorScalarMaterial::nNewMP()),
      _n_new_material_properties_old(
          RankTwoTensorScalarMaterial::nNewMP(/*old=*/true))

{
  unsigned int i;

  // decleare properties
  _new_material_properties.resize(_n_new_material_properties);
  _new_material_properties_old.resize(_n_new_material_properties_old);

  for (i = 0; i < _n_new_material_properties; i++) {
    std::string mp_name = _new_material_property_name;
    if (i == 1) /*eg compute rate flag = true*/
      mp_name += "_rate";
    _new_material_properties[i] = &declareProperty<Real>(mp_name);
  }

  for (i = 0; i < _n_new_material_properties_old; i++) {
    std::string mp_name = _new_material_property_name;
    if (i == 1) /*eg compute rate flag = true*/
      mp_name += "_rate";
    _new_material_properties_old[i] =
        &getMaterialPropertyOldByName<Real>(mp_name);
  }
}

void RankTwoTensorScalarMaterial::initQpStatefulProperties() {
  for (unsigned int i = 0; i < _n_new_material_properties_old; i++)
    (*_new_material_properties[i])[_qp] = 0;
}

void RankTwoTensorScalarMaterial::computeQpProperties() {
  unsigned int qp = check_qp(_qp);

  (*_new_material_properties[0])[_qp] = RankTwoScalarTools::getQuantity(
      _input_rank_two_tensor[qp], _scalar_type, _point1, _point2, _q_point[qp],
      _input_direction);
  if (_compute_rate_flag) {
    if (_t_step > 0)
      (*_new_material_properties[1])[_qp] =
          ((*_new_material_properties[0])[qp] -
           (*_new_material_properties_old[0])[qp]) /
          _dt;
    else
      (*_new_material_properties[1])[_qp] = 0;
  }
}

unsigned int RankTwoTensorScalarMaterial::nNewMP(const bool old /*=fasle*/) {

  unsigned int n_new_mp;
  if (!_stateful_flag && !_compute_rate_flag) {
    if (old)
      n_new_mp = 0;
    else
      n_new_mp = 1;

  } else if (_stateful_flag && !_compute_rate_flag) {
    if (old)
      n_new_mp = 1;
    else
      n_new_mp = 1;

  } else if (!_stateful_flag && _compute_rate_flag) {
    if (old)
      n_new_mp = 1;
    else
      n_new_mp = 2;

  } else if (_stateful_flag && _compute_rate_flag) {
    if (old)
      n_new_mp = 2;
    else
      n_new_mp = 2;

  } else
    mooseError("RankTwoTensorScalarMaterial Error: "
               "somewthing is wrong in checking _stateful_flag && "
               "_compute_rate_flag ");
  return n_new_mp;
}

unsigned int
RankTwoTensorScalarMaterial::check_qp(const unsigned int qp) const {
  unsigned int qp_new = qp;
  if (_has_selected_qp) {
    qp_new = _selected_qp;
    if (_selected_qp >= _q_point.size()) {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      mooseError("RankTwoScalarAux.  selected_qp specified as ", _selected_qp,
                 " but there are only ", _q_point.size(),
                 " quadpoints in the element");
    }
  }
  return qp_new;
}
