//
// Copyright (C) 2013-2018 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//
// Preserve the shared form's existing translation context.
pragma Translator: "ParametricSurvivalAnalysis"

import QtQuick
import QtQuick.Layouts
import JASP.Controls
import JASP

Section
{
	id: modelSection
	property bool mixture: false
	readonly property int modelCount: modelTerms.count
	readonly property int variableCount: modelTerms.variableCount
	readonly property bool selectionAllowsMerging: modelTerms.count == 1 || (modelTerms.count > 1 && interpretModel.value !== "bestAic" && interpretModel.value !== "bestBic")

	title: qsTr("Model")

	CheckBox
	{
		id: advancedSpecificationControl
		name: "advancedSpecification"
		label: qsTr("Advanced specification")
		info: qsTr("Choose a transformation for each main-effect variable. Interactions use the same variable transformations. Transformed factors use their numeric level labels as stress values; factors without a transformation remain categorical.")
	}

	ParametricModelTerms
	{
		name:				"modelTerms"
		id:					modelTerms
		advancedSpecification: advancedSpecificationControl.checked
	}

	DropDown
	{
		id:					interpretModel
		name:				"interpretModel"
		label:				qsTr("Interpret model")
		enabled:			modelTerms.count > 1
		onCountChanged:		if (!(value === "bestAic" || value === "bestBic" || value === "all")) currentIndex = count - 1
		info: mixture ? qsTranslate("ParametricMixtureSurvivalAnalysis", "Select the model to interpret. Defaults to the last specified model. Alternatives are 'All' which produces results for all of the specified models or 'Best' which produces results for the best fitting model based on either AIC or BIC. The selection proceeds within each subgroup from the distribution to the number of components to the model: 'Best' keeps the level of the best fitting model across all remaining distributions, numbers of components, and models, and 'All' selects separately within each of its levels (e.g., all distributions with the best number of components and the best model within each distribution).") : qsTr("Select the model to interpret. Defaults to the last specified model. Alternatives are 'All' which produces results for all of the specified models or 'Best' which produces results for the best fitting model based on either AIC or BIC. If distribution and model selection is specified simultanously, the best model within the best performing distribution is going to be selected. If model selection is specified while all distributions are selected, the best model within each distribution is going to be selected.")
		startValue:			"model1"
		source:
		[
			{
				values: [
					{label: qsTr("All"),		value: "all"},
					{label: qsTr("Best AIC"),	value: "bestAic"},
					{label: qsTr("Best BIC"),	value: "bestBic"}
				]
			},
			{
				values: modelTerms.modelTitles
			}
		]
	}
}
