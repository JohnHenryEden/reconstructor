document.getElementById("confidenceScoreInput").addEventListener("input", function () {
    let value = this.value.trim();

    if (value === "") {
        this.value = ""; // Allow blank (null)
    } else {
        let numValue = parseInt(value, 10);
        if (isNaN(numValue) || numValue < 1 || numValue > 4) {
            this.value = ""; // Reset invalid inputs
        }
    }
});
async function checkAndSaveReaction() {
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    const userId = sessionStorage.getItem('userID');

    if (!reactionId) {
        alert('Reaction not created.');
        return;
    }
    try {
        let response = await fetch(alreadySavedURL, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ user_id: userId, reaction_id: reactionId })
        });

        const data = await response.json();

        if (data.is_reaction_saved) {
            alert('Reaction is already saved.');
        } else {
            var modal = document.getElementById('saveReactionModal');
            modal.style.display = 'block';
            document.getElementById('reactionDescriptionInput').value = sessionStorage.getItem('lastTemplateDescription') || '';
            document.getElementById('modalBackground').style.display = 'block';

            if (!flagsLoaded) {
                loadFlags();
            }
        }
    } catch (error) {
        console.error('Error checking reaction:', error);
    }
}
async function reactionNameExists(shortName, userId) {
    try {
        const response = await fetch(reactionNameExistsURL, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ short_name: shortName, user_id: userId })
        });

        const data = await response.json();
        return data.is_name_saved;
    } catch (error) {
        console.error('Error checking reaction name:', error);
        return false;
    }
}
document.addEventListener('DOMContentLoaded', function () {
    const saveReactionButton = document.getElementById('saveReactionButton');

    document.getElementById('submitSaveReaction').addEventListener('click', async function () {
        const userID = sessionStorage.getItem('userID');
        const urlParams = new URLSearchParams(window.location.search);
        const reactionId = urlParams.get('reaction_id');
        const shortNameInput = document.getElementById('reactionNameInput');
        const shortName = shortNameInput.value;
        const descriptionInput = document.getElementById('reactionDescriptionInput'); // New description input
        const reactionDescription = descriptionInput ? descriptionInput.value : ''; // Get description
        const flagNameElement = document.getElementById('selectedOption');
        let flagName = flagNameElement.textContent.trim();
        const flagIcon = flagNameElement.querySelector('i');
        let flagColor = flagIcon ? flagIcon.style.color : '';
    
        // Convert the color if available
        flagColor = flagColor ? rgbToHex(flagColor) : '';

        // Get Confidence Score
        const confidenceScoreInput = document.getElementById("confidenceScoreInput");
        let confidenceScore = confidenceScoreInput.value.trim() === "" ? null : parseInt(confidenceScoreInput.value);

        // Ensure only valid numbers are sent (1-4) OR null
        if (confidenceScore !== null && (confidenceScore < 1 || confidenceScore > 4)) {
            alert("Confidence Score must be between 1 and 4.");
            return;
        }

        shortNameInput.setCustomValidity(''); // Clear any previous custom validity message
    
        if (userID && reactionId) {
            if (!shortName) {
                alert('Please enter a short name for the reaction.');
                shortNameInput.setCustomValidity('Please enter a short name for the reaction.');
                shortNameInput.reportValidity();
                return; // Prevent form submission
            }
    
            const data = new FormData();
            data.append('userID', userID);
            data.append('reaction_id', reactionId);
            data.append('short_name', shortName);
            data.append('description', reactionDescription); // Include description
            data.append('flag_name', flagName);
            data.append('flag_color', flagColor);
            if (confidenceScore !== null) {
                data.append("confidence_score", confidenceScore); // Include confidence score
            }
            const nameExists = await reactionNameExists(shortName, userID);
            if (nameExists) {
                const userConfirmation = confirm(
                    'Reaction name already exists. Are you sure you want to continue? \nYou will have two reactions with the same name.'
                );
                if (!userConfirmation) {
                    return; // Prevent form submission and keep the modal open
                }
            }            
            fetch(saveReaction, {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken
                },
                body: data,
            })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    alert("Reaction saved successfully!");
                    document.getElementById('saveReactionModal').style.display = 'none';
                    document.getElementById('modalBackground').style.display = 'none';
                } else {
                    alert("Error: " + (data.message || "Failed to save the reaction."));
                }
            })
            .catch(error => console.error('Error:', error));
        } else if (!reactionId) {
            alert("Create the reaction first.");
        } else {
            alert("Please log in.");
        }
    });

    saveReactionButton.addEventListener('click', checkAndSaveReaction);
});
