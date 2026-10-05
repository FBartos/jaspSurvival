pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

BasicThreeButtonTableView
{
	id: observations
	property var predictionPredictors: []
	property bool initialized: false
	property var previousPredictors: []

	name: "customPredictionData"
	modelType: JASP.Simple
	itemType: JASP.String
	initialColumnCount: Math.max(1, predictionPredictors.length)
	initialRowCount: 1
	columnNames: predictionPredictors.length > 0 ? predictionPredictors : [qsTr("Observation")]
	defaultValue: ""
	tableView.parseDefaultValue: false
	tableView.showAddRemoveButtons: false
	buttonsInRow: true
	showButtons: false
	Layout.minimumWidth: observationButtons.implicitWidth
	tableView.anchors.top: observationButtons.bottom
	tableView.anchors.topMargin: jaspTheme.generalAnchorMargin
	preferredHeight: Math.min(300 * jaspTheme.uiScale, tableView.y + tableView.tableHeight)
	buttonAddText: qsTr("Add observation")
	buttonDeleteText: qsTr("Delete observation")
	buttonDeleteEnabled: tableView.rowCount > 1
	onAddClicked: tableView.addRow()
	onDeleteClicked: tableView.removeARow()
	onResetClicked: tableView.reset()
	info: qsTr("Enter one observation per row. Continuous predictors accept numeric expressions. Leave a cell empty to use the mean of that predictor in the fitting data, or the first fitted level for a categorical predictor. Values are entered before any model transformations. Subgroup columns select the corresponding subgroup model.")

	RowLayout
	{
		id: observationButtons
		width: parent.width
		spacing: jaspTheme.columnGroupSpacing
		RoundedButton
		{
			text: observations.buttonAddText
			Layout.fillWidth: true
			Layout.preferredWidth: 1
			Layout.minimumWidth: implicitWidth
			onClicked: { forceActiveFocus(); observations.addClicked(); }
		}
		RoundedButton
		{
			text: observations.buttonDeleteText
			enabled: observations.buttonDeleteEnabled
			Layout.fillWidth: true
			Layout.preferredWidth: 1
			Layout.minimumWidth: implicitWidth
			onClicked: { forceActiveFocus(); observations.deleteClicked(); }
		}
		RoundedButton
		{
			text: observations.buttonResetText
			enabled: observations.buttonResetEnabled
			Layout.fillWidth: true
			Layout.preferredWidth: 1
			Layout.minimumWidth: implicitWidth
			onClicked: { forceActiveFocus(); observations.resetClicked(); }
		}
	}

	function getEditable(columnIndex, rowIndex) { return predictionPredictors.length > 0; }

	// Keep observations attached to their predictor when columns are added or reordered.
	onPredictionPredictorsChanged:
	{
		if (!initialized)
			return;
		let values = [];
		for (let col = 0; col < previousPredictors.length; col++)
		{
			let column = [];
			for (let row = 0; row < tableView.rowCount; row++)
				column.push(tableView.model.data(tableView.model.index(row, col), Qt.DisplayRole));
			values.push(column);
		}
		tableView.columnCount = Math.max(1, predictionPredictors.length);
		for (let col = 0; col < tableView.columnCount; col++)
		{
			let previous = previousPredictors.indexOf(predictionPredictors[col]);
			for (let row = 0; row < tableView.rowCount; row++)
				tableView.itemChanged(col, row, previous < 0 ? "" : String(values[previous][row]), "string");
		}
		previousPredictors = predictionPredictors.slice();
	}
	onTableViewCompleted:
	{
		previousPredictors = predictionPredictors.slice();
		initialized = true;
	}
}
