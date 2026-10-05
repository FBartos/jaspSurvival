pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

Group
{
	id: custom
	property var predictionPredictors: []
	columns: 2
	Layout.fillWidth: true
	preferredWidth: width

	CheckBox
	{
		id: predictions
		name: "customPredictions"
		label: qsTr("Prediction for new observations")
		Layout.columnSpan: 2
		info: qsTr("Predict survival times for specified predictor values using the selected fitted models. Predictions describe the unconditional lifetime distribution, including for observations whose outcomes were missing during fitting.")
	}

	RadioButtonGroup
	{
		id: input
		name: "customPredictionSource"
		title: qsTr("Source")
		Layout.fillWidth: true
		Layout.preferredWidth: (custom.width - custom.columnSpacing) / 2
		Layout.alignment: Qt.AlignTop
		enabled: predictions.checked
		RadioButton
		{
			value: "dataset"
			label: qsTr("Dataset (missing outcomes)")
			checked: true
		}
		RadioButton
		{
			value: "manual"
			label: qsTr("Manual input")
		}
		RadioButton
		{
			value: "csv"
			label: qsTr("CSV file")
		}
		info: qsTr("Dataset predictions include only rows with missing outcome data and at least one specified predictor value. Missing predictors use fitting-data means or the first fitted categorical level. The input table and CSV file contain predictor values on their original scale. CSV column names must match the dataset predictor names; subgroup variables must also be supplied when subgroup models are fitted.")
	}

	Group
	{
		title: qsTr("Predict")
		Layout.fillWidth: true
		Layout.preferredWidth: (custom.width - custom.columnSpacing) / 2
		Layout.alignment: Qt.AlignTop
		enabled: predictions.checked
		CheckBox
		{
			name: "customPredictionMean"
			label: qsTr("Mean")
			checked: true
			info: qsTr("Predict the mean survival time for each observation.")
		}
		CheckBox
		{
			name: "customPredictionPercentiles"
			label: qsTr("Percentiles")
			childrenOnSameRow: true
			info: qsTr("Predict the survival times below which the specified percentages of events occur.")
			FormulaField
			{
				name: "customPredictionPercentileValues"
				defaultValue: "25, 50, 75"
				inputType: "string"
				parseDefaultValue: false
				fieldWidth: 150 * jaspTheme.uiScale
				info: qsTr("Specify percentages from 0 to 100, separated by commas or semicolons. For example, 50 predicts the median survival time.")
			}
		}
		CheckBox
		{
			name: "customPredictionConfidenceInterval"
			label: qsTr("Confidence intervals")
			checked: true
			childrenOnSameRow: true
			info: qsTr("Include confidence intervals for the predicted mean and percentiles. These describe uncertainty in the fitted quantities.")
			FormulaField
			{
				name: "customPredictionConfidenceLevels"
				defaultValue: "95"
				inputType: "string"
				parseDefaultValue: false
				afterLabel: "%"
				// Match CIField's standard width; retain entry of multiple levels.
				fieldWidth: jaspTheme.font.pixelSize * 4
				info: qsTr("Specify one or more confidence levels between 0 and 100, separated by commas or semicolons, such as 90, 95, 99.")
			}
		}
	}
	ParametricPredictionData
	{
		Layout.columnSpan: 2
		Layout.fillWidth: true
		predictionPredictors: custom.predictionPredictors
		visible: input.value === "manual"
		enabled: predictions.checked
	}
	FileSelector
	{
		Layout.columnSpan: 2
		name: "customPredictionFile"
		label: qsTr("CSV file")
		save: false
		filter: "*.csv"
		fieldWidth: 250 * jaspTheme.uiScale
		visible: input.value === "csv"
		enabled: predictions.checked
		info: qsTr("Select a comma-separated CSV file with a header row and one observation per row. Predictor names and categorical levels must match the fitting dataset. Extra columns are ignored.")
	}
	Group
	{
		Layout.columnSpan: 2
		title: qsTr("Output")
		enabled: predictions.checked
		CheckBox
		{
			name: "customPredictionTable"
			label: qsTr("Table")
			checked: true
			info: qsTr("Display predictions and the predictor values for each observation.")
		}
		CheckBox
		{
			name: "customPredictionExport"
			label: qsTr("Export predictions to dataset")
			checked: true
			enabled: input.value === "dataset"
			info: qsTr("Add predictions only for rows with missing outcomes and at least one specified predictor value. Missing predictors use fitting-data defaults. Other rows and unsupported predictor values remain missing in the prediction columns.")
			TextField
			{
				name: "customPredictionColumnPrefix"
				label: qsTr("Column prefix")
				defaultValue: "Prediction"
				fieldWidth: 150 * jaspTheme.uiScale
				info: qsTr("Prefix for exported prediction columns. Model identifiers are appended when multiple models are selected.")
			}
		}
	}
}
