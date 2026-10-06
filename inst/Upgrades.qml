import QtQuick
import JASP.Module

Upgrades
{
	Upgrade
	{
		functionName: "ParametricSurvivalAnalysis"
		fromVersion: "0.96.7"
		toVersion: "0.97.0"
		msg: qsTr("Parametric survival models now always include an intercept. Results may change if this saved analysis previously excluded the intercept.")

		ChangeRemove
		{
			name: "includeIntercept"
			condition: function(options) { return options["includeIntercept"] !== undefined; }
		}

		ChangeJS
		{
			name: "subgroup"
			jsFunction: function(options)
			{
				let subgroup = options["subgroup"];
				if (typeof subgroup === "string")
					return subgroup.length > 0 ? [subgroup] : [];
				if (subgroup && typeof subgroup.value === "string")
					subgroup.value = subgroup.value.length > 0 ? [subgroup.value] : [];
				return subgroup;
			}
		}
	}

	Upgrade
	{
		functionName: "ParametricSurvivalAnalysis"
		fromVersion: "0.97.0"
		toVersion: "0.97.1"

		ChangeSetValue
		{
			name: "advancedSpecification"
			condition: function(options) { return options["advancedSpecification"] === undefined; }
			jsonValue: false
		}
		ChangeSetValue
		{
			name: "restrictParameters"
			condition: function(options) { return options["restrictParameters"] === undefined; }
			jsonValue: false
		}
		ChangeSetValue
		{
			name: "fixedParameters"
			condition: function(options) { return options["fixedParameters"] === undefined; }
			jsonValue: []
		}
		ChangeSetValue
		{
			name: "customPredictions"
			condition: function(options) { return options["customPredictions"] === undefined; }
			jsonValue: false
		}
	}

	Upgrade
	{
		functionName: "ParametricMixtureSurvivalAnalysis"
		fromVersion: "0.97.0"
		toVersion: "0.97.1"

		ChangeSetValue
		{
			name: "advancedSpecification"
			condition: function(options) { return options["advancedSpecification"] === undefined; }
			jsonValue: false
		}
		ChangeSetValue
		{
			name: "restrictParameters"
			condition: function(options) { return options["restrictParameters"] === undefined; }
			jsonValue: false
		}
		ChangeSetValue
		{
			name: "fixedParameters"
			condition: function(options) { return options["fixedParameters"] === undefined; }
			jsonValue: []
		}
		ChangeSetValue
		{
			name: "customPredictions"
			condition: function(options) { return options["customPredictions"] === undefined; }
			jsonValue: false
		}
	}
}
