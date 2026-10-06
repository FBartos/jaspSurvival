pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Controls as QtControls
import QtQuick.Layouts
import JASP.Controls
import JASP

// Rich labels keep mathematical subscripts readable in both the selection and menu.
DropDown
{
	id: transformation
	property bool allowScale: false
	Layout.alignment: Qt.AlignTop
	name: "transformation"
	startValue: "none"
	fieldWidth: 150 * jaspTheme.uiScale
	useExternalBorder: true
	control.contentItem: Text
	{
		text: transformation.currentLabel
		textFormat: Text.RichText
		font: transformation.control.font
		color: enabled ? jaspTheme.textEnabled : jaspTheme.textDisabled
		verticalAlignment: Text.AlignVCenter
		leftPadding: 4 * jaspTheme.uiScale
		clip: true
	}
	control.delegate: QtControls.ItemDelegate
	{
		id: transformationItem
		required property int index
		required property var model
		width: transformation ? transformation.control.popup.width : 0
		height: jaspTheme.comboBoxHeight
		readonly property bool selected: transformation && transformation.currentIndex === index
		highlighted: transformation && transformation.control.highlightedIndex === index
		contentItem: Text
		{
			text: model.name
			textFormat: Text.RichText
			font: jaspTheme.font
			color: transformationItem.selected ? jaspTheme.white : jaspTheme.textEnabled
			verticalAlignment: Text.AlignVCenter
		}
		background: Rectangle
		{
			color: transformationItem.selected ? jaspTheme.itemSelectedColor :
				(transformationItem.hovered || transformationItem.highlighted ? jaspTheme.itemHoverColor : "transparent")
		}
	}
	values:
	{
		let transformations = [
			{label: qsTr("None"), value: "none"},
			{label: qsTr("Power: ln(x)"), value: "power"},
			{label: qsTr("Power (y): log<sub>y</sub>(x)"), value: "customPower"},
			{label: qsTr("Inverse power: -ln(x)"), value: "inversePower"},
			{label: qsTr("Inverse power (y): -log<sub>y</sub>(x)"), value: "customInversePower"},
			{label: qsTr("Arrhenius: 1/T"), value: "arrhenius"},
			{label: qsTr("Reciprocal: 1/x"), value: "reciprocal"},
			{label: qsTr("Square: x²"), value: "square"},
			{label: qsTr("Exponent (y): xʸ"), value: "customExponent"},
			{label: qsTr("Square root: √x"), value: "squareRoot"},
			{label: qsTr("Polynomial"), value: "polynomial"}
		];
		if (allowScale)
			transformations.splice(1, 0, {label: qsTr("Scale: (x - mean(x)) / sd(x)"), value: "scale"});
		return transformations;
	}
	info: qsTr("Transform predictors before applying the distribution's model link. None keeps the previous model and preserves variable types. Scale centers covariates at their mean and divides by their sample standard deviation, using the data fitted by each model. Power and inverse power use natural logarithms; Power (y) and Inverse power (y) use the selected base. Square, square root, and reciprocal are fixed powers; Exponent (y) uses the selected exponent. Arrhenius uses reciprocal absolute temperature. Polynomial coefficients refer to orthogonal polynomial terms.")
}
