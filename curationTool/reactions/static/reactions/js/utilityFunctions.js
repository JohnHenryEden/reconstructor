// Generates a list of distinct colors in HSL format, avoiding yellow hues, based on the specified count.
function getDistinctColors(count) {
    var colors = [];
    var hueStep = 360 / count;
    var lightness = 40; // Adjust for darker colors, but not too dark

    for (var i = 0; i < 360; i += hueStep) {
        // Skip yellow and near-yellow hues
        if (i >= 50 && i <= 70) {
            continue;
        }

        colors.push('hsl(' + i + ', 100%, ' + lightness + '%)');
    }
    return colors;
}

window.ReactionUtils = window.ReactionUtils || {};

window.ReactionUtils.fetchTemplates = async function (templateList, templatesFetched) {
    if (templatesFetched.value) return;
    try {
        const userID = sessionStorage.getItem('userID');
        console.log(userID);
        const response = await fetch('/list_templates/', {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ userID: userID }),
        });

        if (response.ok) {
            const data = await response.json();
            templateList.length = 0; // Clear existing list
            templateList.push(...data.templates); // Add new templates
            templatesFetched.value = true;
        } else {
            console.error('Failed to fetch templates');
        }
    } catch (error) {
        console.error('Error:', error);
    }
};
