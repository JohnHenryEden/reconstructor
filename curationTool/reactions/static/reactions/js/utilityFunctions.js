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

// Function to show a toast notification
function showToast(message) {
    // Create the toast container
    const toast = document.createElement('div');
    toast.textContent = message;
    toast.style.position = 'fixed';
    toast.style.bottom = '20px';
    toast.style.right = '20px';
    toast.style.backgroundColor = '#4caf50'; // Success green
    toast.style.color = '#fff';
    toast.style.padding = '10px 20px';
    toast.style.borderRadius = '8px';
    toast.style.boxShadow = '0 4px 8px rgba(0, 0, 0, 0.2)';
    toast.style.fontSize = '14px';
    toast.style.zIndex = '1000';
    toast.style.opacity = '0';
    toast.style.transition = 'opacity 0.5s ease';

    document.body.appendChild(toast);

    // Fade in the toast
    setTimeout(() => {
        toast.style.opacity = '1';
    }, 10);

    // Fade out and remove the toast after 3 seconds
    setTimeout(() => {
        toast.style.opacity = '0';
        setTimeout(() => toast.remove(), 500); // Remove after fade-out
    }, 3000);
}
