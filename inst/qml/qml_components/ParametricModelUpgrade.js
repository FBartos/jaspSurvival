// Upgrade the old FactorsForm option to bound per-term rows.
function modelTerms(options)
{
	return options["modelTerms"].map(function(model)
	{
		let components = model.components;
		let typed = components && components.value !== undefined;
		let values = typed ? components.value : components;
		let types = typed ? components.types : [];
		let rows = values.map(function(term)
		{
			return {
				components: Array.isArray(term) ? term : [term],
				transformation: "none",
				temperatureUnit: "kelvin"
			};
		});
		model.components = typed ? {value: rows, types: types} : rows;
		return model;
	});
}
