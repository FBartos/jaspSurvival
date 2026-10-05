// Bound model lists retain the per-variable controls in saved analyses.
pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

ComponentsList
{
	id: models
	property bool advancedSpecification: false
	property int variableCount: 0
	property var modelTitles: []
	optionKey: "name"
	newItemValue: "model1"
	minimumItems: 1
	duplicateWhenAdding: true
	addBorder: false
	addTooltip: qsTr("Add a model")
	removeTooltip: qsTr("Remove this model")
	info: qsTr("Specify candidate models. Each variable has one transformation within a model, shared by its main effect and interactions. Add a main effect to specify its transformation. Models can use different transformations and need not be nested.")

	function updateModels()
	{
		let titles = [];
		let variables = 0;
		for (let key of columnsNames)
		{
			let title = getRowControl(key, "title");
			let terms = getRowControl(key, "components");
			if (title) titles.push({label: title.value, value: key});
			if (terms) variables += terms.count;
		}
		modelTitles = titles;
		variableCount = variables;
	}
	onColumnsNamesChanged: updateModels()
	onBoundValueChanged: updateModels()

	rowComponent: Group
	{
		TextField
		{
			name: "title"
			defaultValue: qsTr("Model %1").arg(rowIndex + 1)
			info: qsTr("Name of this candidate model.")
			onValueChanged: models.updateModels()
		}

		VariablesForm
		{
			id: modelVariables
			preferredHeight: 180 * jaspTheme.uiScale
			AvailableVariablesList
			{
				id: availableModelTerms
				name: "availableModelTerms"
				width: models.advancedSpecification ? 200 * jaspTheme.uiScale : modelVariables.listWidth
				source: ['covariates', 'factors']
				keepVariablesWhenMoved: true
			}
			AssignedVariablesList
			{
				id: terms
				name: "components"
				width: models.advancedSpecification ? modelVariables.width - availableModelTerms.width - 80 * jaspTheme.uiScale : modelVariables.listWidth
				title: qsTr("Model terms")
				listViewType: JASP.Interaction
				allowedColumns: []
				addAvailableVariablesToAssigned: true
				onCountChanged: models.updateModels()
				rowComponentTitle: models.advancedSpecification ? qsTr("Relationship / transformation") : ""
				rowComponent: RowLayout
				{
					visible: models.advancedSpecification
					// Interaction rows use the choices on their constituent main effects.
					readonly property bool mainEffect: rowType !== "unknown"
					DropDown
					{
						name: "transformation"
						id: transformation
						visible: parent.mainEffect
						startValue: "none"
						values: [
							{label: qsTr("None"), value: "none"},
							{label: qsTr("Exponential: x"), value: "exponential"},
							{label: qsTr("Power: ln(x)"), value: "power"},
							{label: qsTr("Inverse power: -ln(x)"), value: "inversePower"},
							{label: qsTr("Arrhenius: 1/T"), value: "arrhenius"},
							{label: qsTr("Reciprocal: 1/x"), value: "reciprocal"},
							{label: qsTr("Square: x²"), value: "square"},
							{label: qsTr("Square root: √x"), value: "squareRoot"}
						]
						info: qsTr("Transform the numeric stress values before fitting. None preserves the variable's original type. Exponential uses numeric values; power and inverse power use natural logarithms. Arrhenius uses reciprocal absolute temperature. The relationship coefficient is estimated freely.")
					}
					DropDown
					{
						name: "temperatureUnit"
						visible: parent.mainEffect && transformation.value === "arrhenius"
						startValue: "kelvin"
						values: [
							{label: qsTr("Kelvin"), value: "kelvin"},
							{label: qsTr("Celsius"), value: "celsius"}
						]
						info: qsTr("Temperature units in the data. Celsius is converted to kelvin before taking the reciprocal.")
					}
					Label
					{
						visible: !parent.mainEffect
						text: qsTr("Uses main-effect transformations")
					}
				}
			}
		}
		Component.onCompleted: models.updateModels()
	}
}
