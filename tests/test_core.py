import openmc_depletion_plotter


def test_read_depletion_results():
    results = openmc_depletion_plotter.read_depletion_results('tests/depletion_results.h5')
    assert len(results) == 25

    # this is the unirradiated material so it should return just Fe56
    condition = results["timestep"] == 0
    filter_data = results[condition][["material_id", "timestep", "nuclide", "atoms"]]
    assert len(filter_data)  == 1
    assert filter_data[0].timestep == 0
    assert filter_data[0].nuclide == 'Fe56'
    assert filter_data[0].atoms == 4.7307702602013954e29
    assert filter_data[0].material_id == 1

    # filter should return lots of results as this is just after irradiation
    condition = results["timestep"] == 3
    filter_data = results[condition][["material_id", "timestep", "nuclide", "atoms"]]
    assert len(filter_data) > 1

    # filter should return no matching results as there are only 3 timesteps
    condition = results["timestep"] == 4
    filter_data = results[condition][["material_id", "timestep", "nuclide", "atoms"]]
    assert len(filter_data) == 0