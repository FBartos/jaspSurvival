pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

ComponentsList
{
	id: modifiers
	property var models
	property bool mixture: false
	property int componentCount: 1
	depends: models
	source: [{values: models.factorsTitles}]
	optionKey: "model"
	addItemManually: false
	addBorder: false
	rowSpacing: 10 * jaspTheme.uiScale
	info: qsTr("Specify regression coefficient restrictions and term transformations separately for each candidate model. Interactions use the transformations of their main effects.")

	rowComponent: ComponentsList
	{
		property bool mixture: listView.mixture
		property int componentCount: listView.componentCount
		implicitWidth: listView.availableWidth
		name: "predictors"
		source: [{name: rowValue, use: "noInteraction"}]
		addItemManually: false
		addBorder: true
		itemRectangle.color: jaspTheme.controlBackgroundColor
		optionKey: "variable"
		preferredWidth: implicitWidth
		rowSpacing: 4 * jaspTheme.uiScale
		headerLabels: [{variableColumn: rowLabel}].concat(
			mixture ? [{componentColumn: qsTr("Component")}, {componentRestrictions: qsTr("Parameter restriction")}] : [{parameterRestriction: qsTr("Parameter restriction")}],
			[{transformation: qsTr("Term transformation")}])

		rowComponent: RowLayout
		{
			id: predictorRow
			property bool mixture: listView.mixture
			property int componentCount: listView.componentCount
			spacing: 8 * jaspTheme.uiScale
			Group
			{
				name: "variableColumn"
				Layout.preferredWidth: 150 * jaspTheme.uiScale
				Layout.alignment: predictorRow.mixture ? Qt.AlignTop : Qt.AlignVCenter
				Label
				{
					text: rowLabel
					height: predictorRow.mixture ? componentRestrictions.rowHeight : implicitHeight
					verticalAlignment: Text.AlignVCenter
				}
			}
			Group
			{
				name: "componentColumn"
				visible: predictorRow.mixture
				Layout.preferredWidth: visible ? 60 * jaspTheme.uiScale : 0
				Layout.alignment: predictorRow.mixture ? Qt.AlignTop : Qt.AlignVCenter
				Column
				{
					spacing: componentRestrictions.rowSpacing
					Repeater
					{
						model: predictorRow.componentCount
						Label
						{
							text: (index + 1).toString()
							height: componentRestrictions.rowHeight
							verticalAlignment: Text.AlignVCenter
						}
					}
				}
			}
			ParametricFixedParameter
			{
				name: "parameterRestriction"
				visible: !predictorRow.mixture
				positive: false
				fieldWidth: 125 * jaspTheme.uiScale
				info: qsTr("Fix main-effect regression coefficients to numeric values or expressions; leave empty to estimate them. For factors, provide one value for every coefficient or a comma- or semicolon-separated vector matching the coefficients shown in the results table (one fewer than the number of levels). These restrictions apply to this model in every selected distribution and mixture component. Transformed factors are treated as numeric predictors.")
			}
			ParametricRegressionComponentRestrictions
			{
				id: componentRestrictions
				name: "componentRestrictions"
				componentCount: predictorRow.componentCount
				visible: predictorRow.mixture
			}
			ParametricTermTransformation
			{
				id: transformation
				allowScale: rowType === "scale"
			}
			DoubleField
			{
				name: "customLogBase"
				label: qsTr("Base")
				visible: transformation.value === "customPower" || transformation.value === "customInversePower"
				enabled: transformation.value === "customPower" || transformation.value === "customInversePower"
				defaultValue: 10
				min: 0
				inclusive: JASP.MaxOnly
				decimals: 6
				fieldWidth: 40 * jaspTheme.uiScale
				useExternalBorder: true
				info: qsTr("Base y of the logarithm. The default is 10. The base must be positive and different from 1.")
			}
			DoubleField
			{
				name: "powerExponent"
				label: qsTr("Exponent")
				visible: transformation.value === "customExponent"
				enabled: transformation.value === "customExponent"
				defaultValue: 2
				negativeValues: true
				decimals: 6
				fieldWidth: 40 * jaspTheme.uiScale
				useExternalBorder: true
				info: qsTr("Exponent y in x raised to y. The default is 2. Negative and fractional exponents are allowed. The exponent must be different from 0.")
			}
			IntegerField
			{
				name: "polynomialDegree"
				label: qsTr("Degree")
				visible: transformation.value === "polynomial"
				enabled: transformation.value === "polynomial"
				defaultValue: 2
				min: 1
				max: 3
				fieldWidth: 40 * jaspTheme.uiScale
				useExternalBorder: true
				info: qsTr("Polynomial degree: 1 fits a linear curve, 2 a quadratic curve, and 3 a cubic curve. Orthogonal polynomial terms improve numerical stability. Interactions include all selected polynomial terms.")
			}
			DropDown
			{
				name: "temperatureUnit"
				visible: transformation.value === "arrhenius"
				startValue: "kelvin"
				fieldWidth: 90 * jaspTheme.uiScale
				useExternalBorder: true
				values: [
					{label: qsTr("Kelvin"), value: "kelvin"},
					{label: qsTr("Celsius"), value: "celsius"}
				]
				info: qsTr("Temperature units in the data. Celsius is converted to kelvin before taking the reciprocal.")
			}
		}
	}
}
