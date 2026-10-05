pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

ComponentsList
{
	id: restrictions
	property var selector
	name: "fixedParameters"
	optionKey: "restriction"
	addItemManually: false
	addBorder: true
	itemRectangle.color: jaspTheme.controlBackgroundColor
	visible: selector.parametersRestricted
	enabled: selector.parametersRestricted
	Layout.columnSpan: form.columns
	Layout.fillWidth: true
	rowSpacing: 4 * jaspTheme.uiScale
	headerLabels: [{distributionColumn: qsTr("Distribution")}].concat(
		selector.mixture ? [{componentColumn: qsTr("Component")}] : [],
		[{parametersColumn: qsTr("Parameters")}])
	readonly property var restrictionRows: {
		let families = ["all", "bestAic", "bestBic"].indexOf(selector.value) >= 0 ?
			selector.selectedFamilies : [{label: selector.currentLabel, value: selector.value}];
		let rows = [];
		for (let family of families)
			for (let k = 1; k <= selector.componentCount; k++)
				rows.push({label: family.label, value: family.value + ":" + k});
		return rows;
	}
	source: [{values: restrictionRows}]

	rowComponent: RowLayout
	{
		readonly property var specification: rowValue.split(":")
		spacing: 8 * jaspTheme.uiScale

		Group
		{
			name: "distributionColumn"
			Layout.preferredWidth: 190 * jaspTheme.uiScale
			Label { text: specification[1] === "1" ? rowLabel : "" }
		}
		Group
		{
			name: "componentColumn"
			visible: restrictions.selector.mixture
			Layout.preferredWidth: visible ? 60 * jaspTheme.uiScale : 0
			Label
			{
				text: specification[1]
				visible: restrictions.selector.mixture
			}
		}
		Group
		{
			name: "parametersColumn"
			ParametricFixedParameterFields { distribution: specification[0] }
		}
	}
}
