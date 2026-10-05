pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

ComponentsListBase
{
	id: restrictions
	property var selector
	property var restrictionKeys: []
	property int sourceGeneration: 0
	property real rowSpacing: 4 * jaspTheme.uiScale
	property bool addBorder: true
	property alias itemRectangle: frame
	readonly property real distributionWidth: 190 * jaspTheme.uiScale
	readonly property real componentWidth: 60 * jaspTheme.uiScale
	readonly property var activeFamilies: ["all", "bestAic", "bestBic"].indexOf(selector.value) >= 0 ?
		selector.selectedFamilies.map(function(family) { return family.value; }) : [selector.value]
	name: "fixedParameters"
	optionKey: "restriction"
	addItemManually: false
	visible: selector.parametersRestricted
	enabled: selector.parametersRestricted
	Layout.columnSpan: form.columns
	Layout.fillWidth: true
	implicitHeight: rows.height + 2 * jaspTheme.contentMargin
	preferredWidth: parent.width
	preferredHeight: implicitHeight
	background: frame

	// Stable slots avoid restoring old saved values when component counts or families change.
	source: [{values: restrictionKeys.map(function(key) {
		return {value: key, label: sourceGeneration + ":" + key};
	})}]
	Component.onCompleted: {
		let keys = [];
		for (let family of selector.familyChoices)
			for (let k = 1; k <= 4; k++) keys.push(family.value + ":" + k);
		restrictionKeys = keys;
	}
	// Older saved options contain only the active slots; populate missing slots once after loading.
	onCountChanged: if (count < restrictionKeys.length) expandSlots.restart()
	Timer
	{
		id: expandSlots
		interval: 0
		onTriggered: if (restrictions.count < restrictions.restrictionKeys.length) restrictions.sourceGeneration++
	}

	Rectangle
	{
		id: frame
		anchors.fill: parent
		color: jaspTheme.controlBackgroundColor
		border.width: restrictions.addBorder ? 1 : 0
		border.color: jaspTheme.borderColor
		radius: jaspTheme.borderRadius
	}
	Column
	{
		id: rows
		width: parent.width - 2 * jaspTheme.contentMargin
		x: jaspTheme.contentMargin
		y: jaspTheme.contentMargin
		spacing: restrictions.rowSpacing
		RowLayout
		{
			spacing: 8 * jaspTheme.uiScale
			Label { text: qsTr("Distribution"); Layout.preferredWidth: restrictions.distributionWidth }
			Label { text: qsTr("Component"); visible: restrictions.selector.mixture; Layout.preferredWidth: restrictions.componentWidth }
			Label { text: qsTr("Parameters") }
		}
		Repeater
		{
			model: restrictions.model
			delegate: FocusScope
			{
				property var rowComponentItem: model.rowComponent
				readonly property var specification: String(model.value).split(":")
				visible: restrictions.activeFamilies.indexOf(specification[0]) >= 0 && Number(specification[1]) <= restrictions.selector.componentCount
				width: rowComponentItem ? rowComponentItem.width : 0
				height: rowComponentItem ? rowComponentItem.height : 0
				Component.onCompleted: if (rowComponentItem) rowComponentItem.parent = this
			}
		}
	}
	rowComponent: RowLayout
	{
		id: restrictionRow
		readonly property var specification: rowValue.split(":")
		readonly property bool mixture: listView.selector.mixture
		spacing: 8 * jaspTheme.uiScale
		Group
		{
			name: "distributionColumn"
			Layout.preferredWidth: listView.distributionWidth
			Label { text: restrictionRow.specification[1] === "1" ? listView.selector.familyChoices.filter(function(family) { return family.value === restrictionRow.specification[0]; })[0].label : "" }
		}
		Group
		{
			name: "componentColumn"
			visible: restrictionRow.mixture
			Layout.preferredWidth: visible ? listView.componentWidth : 0
			Label { text: restrictionRow.specification[1]; visible: restrictionRow.mixture }
		}
		Group
		{
			name: "parametersColumn"
			ParametricFixedParameterFields { distribution: restrictionRow.specification[0] }
		}
	}
}
