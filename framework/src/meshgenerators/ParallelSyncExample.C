//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParallelSyncExample.h"
#include "InputParameters.h"
#include "libmesh/parallel_sync.h"

registerMooseObject("MooseApp", ParallelSyncExample);

defineLegacyParams(ParallelSyncExample);

InputParameters
ParallelSyncExample::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addClassDescription("paralle sync example");

  return params;
}

ParallelSyncExample::ParallelSyncExample(const InputParameters & parameters)
  : MeshGenerator(parameters), _input(getMesh("input"))
{
}

std::unique_ptr<MeshBase>
ParallelSyncExample::generate()
{

  if (_input->is_replicated())
    mooseError("ParallelSyncExample should only be used in distributed mode");

  std::unique_ptr<MeshBase> mesh = std::move(_input);
  // generate node to element map
  std::map<dof_id_type, std::vector<dof_id_type>> node_to_elem_map;
  for (const auto & elem : mesh->active_element_ptr_range())
    if (!elem->is_remote())
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
        node_to_elem_map[elem->node_id(n)].push_back(elem->id());

  std::vector<dof_id_type> node_multiplicity_local = findNodeMultiplicity(node_to_elem_map, *mesh);
  getMyDataParallel(node_multiplicity_local, *mesh);

  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}

void
ParallelSyncExample::getMyDataParallel(const std::vector<dof_id_type> & node_multiplicity_local,
                                       const MeshBase & mesh)
{

  // exchange this info with other processors, building _nodes_to_send at the same time
  std::map<processor_id_type, std::vector<dof_id_type>> queries, received_data;

  queries.clear();
  received_data.clear();

  // set query to local node multiplicity
  queries[mesh.comm().rank()] = node_multiplicity_local;

  // compose replies from each cpu
  auto gather_data = [](processor_id_type pid,
                        const std::vector<dof_id_type> & node_multiplicity_local,
                        std::vector<dof_id_type> & response) {
    response = node_multiplicity_local;
  };

  // put togherter all teh responses
  auto act_on_data = [&received_data](processor_id_type pid,
                                      const std::vector<dof_id_type> & node_multiplicity_local,
                                      const std::vector<dof_id_type> & response) {
    auto & vec = received_data[pid];
    vec.insert(vec.end(), response.begin(), response.end());
  };

  dof_id_type * ex = nullptr;
  Parallel::pull_parallel_vector_data(mesh.comm(), queries, gather_data, act_on_data, ex);

  // print some output
  for (auto const & x : received_data)
    std::cout << " Rank " << mesh.comm().rank() << "  received data from Rank " << x.first
              << "!!!!!!!!!!!!!!!!!!! " << std::endl;
}

std::vector<dof_id_type>
ParallelSyncExample::findNodeMultiplicity(
    const std::map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map, const MeshBase & mesh)
{
  std::vector<dof_id_type> node_multiplicity_local(0);
  for (auto node_it = node_to_elem_map.begin(); node_it != node_to_elem_map.end(); ++node_it)
  {
    const dof_id_type current_node_id = node_it->first;
    const Node * current_node = mesh.node_ptr(current_node_id);

    if (current_node != nullptr && current_node)
    {
      // find node multiplicity
      std::set<subdomain_id_type> connected_blocks;
      for (auto elem_id = node_it->second.begin(); elem_id != node_it->second.end(); elem_id++)
      {

        const Elem * current_elem = mesh.elem_ptr(*elem_id);
        if (current_elem->processor_id() == comm().rank())
          connected_blocks.insert(current_elem->subdomain_id());
      }
      unsigned int node_multiplicity = connected_blocks.size();
      if (node_multiplicity > 1)
        for (auto block_id : connected_blocks)
        {
          node_multiplicity_local.push_back(current_node_id);
          node_multiplicity_local.push_back(block_id);
        }
    }
  }
  return node_multiplicity_local;
}
