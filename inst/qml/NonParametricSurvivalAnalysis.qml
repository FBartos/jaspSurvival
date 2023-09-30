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
import QtQuick			2.8
import QtQuick.Layouts	1.3
import JASP.Controls	1.0
import JASP.Widgets		1.0
import JASP				1.0

Form
{
	VariablesForm
	{
		AvailableVariablesList
		{
			name: "allVariablesList"
		}

		AssignedVariablesList
		{
			name:				"timeToEvent"
			title:				qsTr("Time to Event")
			suggestedColumns:	["scale"]
			singleVariable:		true
		}

		AssignedVariablesList
		{
			name:				"eventStatus"
			title:				qsTr("Event Status")
			suggestedColumns:	["nominal"]
			singleVariable:		true
		}

		DropDown
		{
			name:				"eventIndicator"
			label:				qsTr("Event Indicator")
			source:				[{name: "eventStatus", use: "levels"}]
		}

		AssignedVariablesList
		{
			name:			 	"factors"
			title:			 	qsTr("Factors")
			allowedColumns:		["ordinal", "nominal", "nominalText"]
		}
	}

	CheckBox
	{
		name:	"lifeTable"
		label:	qsTr("Life table")

		CheckBox
		{
			name:		"lifeTableRoundSteps"
			label:		qsTr("Round steps")
			checked:	true
		}

		DropDown
		{
			name:		"lifeTableStepsType"
			id:			lifeTableStepsType
			label:		qsTr("Steps type")
			values:
			[
				{ label: qsTr("Quantilies"),	value: "quantiles"},
				{ label: qsTr("Fixed size"),	value: "fixedSize"}
			]
		}

		IntegerField
		{
			name:			"lifeTableStepsNumber"
			label:			qsTr("Number")
			defaultValue:	10
			visible:		lifeTableStepsType.value == "quantiles"
		}

		DoubleField
		{
			name:			"lifeTableStepsFrom"
			id:				lifeTableStepsFrom
			label:			qsTr("From")
			defaultValue:	0
			max:			lifeTableStepsTo.value
			visible:		lifeTableStepsType.value == "fixedSize"
		}

		DoubleField
		{
			name:			"lifeTableStepsSize"
			label:			qsTr("Size")
			defaultValue:	1
			// max:			lifeTableStepsTo.value // TODO: enable once max is data dependent
			visible:		lifeTableStepsType.value == "fixedSize"
		}

		DoubleField
		{
			name:			"lifeTableStepsTo"
			id:				lifeTableStepsTo
			label:			qsTr("To")
			defaultValue:	0
			min:			lifeTableStepsFrom.value
			visible:		lifeTableStepsType.value == "fixedSize"
		}
	}
}
