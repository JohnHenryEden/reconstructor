document.getElementById('submitBtn-form').addEventListener('click', function(event) {
    event.preventDefault(); // Prevent the default form submission

    // Trigger the form's submit event
    document.getElementById('reactionForm').requestSubmit();
});

function hidemodal(){

    document.getElementById('error-modal').style.display = 'none';

}

document.getElementById('close-button').onclick = function() {
    document.getElementById('error-modal').style.display = 'none';
};

function showIdenticalReactionModal(matches, submitBtn) {
    return new Promise((resolve) => {
        const modal = $('#identicalReactionModal');

        // 1. Build a dynamic message
        let messageText = '';
        if (matches.length > 1) {
            messageText += `Multiple identical reactions have already been saved in your reactions: <br><br>`;
        } else {
            messageText += `An identical reaction has already been saved in your reactions: <br><br>`;
        }

        // 2. Generate a "View" button for each match
        matches.forEach(({ reaction_id, reaction_name }) => {
            messageText += `
                <button 
                  class="ui button dynamic-view-btn" 
                  data-reaction-id="${reaction_id}" 
                  style="margin-bottom: 0.5em; margin-right: 0.5em;"
                >
                  View "${reaction_name}"
                </button>
            `;
        });

        // 3. Prompt user they can create new
        messageText += `<br><br>Click "Create" to proceed with creating a new reaction.`;

        // Insert the HTML into the modal
        document.getElementById('identicalReactionMessage').innerHTML = messageText;

        // 4. Attach handlers for each generated "View" button
        setTimeout(() => {
            const viewButtons = document.querySelectorAll('.dynamic-view-btn');
            viewButtons.forEach(btn => {
                btn.addEventListener('click', function() {
                    const rxnId = this.getAttribute('data-reaction-id');
                    resolve(`view:${rxnId}`);
                    modal.modal('hide');
                });
            });
        }, 0);

        // 5. Keep the existing "Create" button
        document.getElementById('createReactionButton').onclick = function () {
            resolve('create');
            modal.modal('hide');
        };

        // Re-enable submit button if the user closes modal
        modal.modal({
            onHide: function() {
                submitBtn.disabled = false;
            }
        });

        // Show the modal
        modal.modal('show');
    });
}


async function checkIdenticalReaction(formData, loadingIndicator, submitBtn) {
    try {
        const response = await fetch(identicalReactionUrl, {
            method: 'POST',
            body: formData,
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            }
        });

        const data = await response.json();

        if (data.status === 'error') {
            showErrorModal(data.message);
            window.scrollTo(0, 0);
            submitBtn.disabled = false;
            loadingIndicator.style.display = 'none';
            return false;
        }

        if (data.exists) {
            loadingIndicator.style.display = 'none';
            const userDecision = await showIdenticalReactionModal(data.matches, submitBtn);
            if (userDecision.startsWith('view:')) {
                const rxnId = userDecision.split(':')[1];
                window.location.href = `/?reaction_id=${rxnId}`;
                return false;
            }
        }
        loadingIndicator.style.display = 'flex';
        return true; // Proceed with reaction creation
    } catch (error) {
        console.error('Error checking for identical reaction:', error);
        return false;
    }
}

document.getElementById('reactionForm').addEventListener('submit', async function(e) {
    var loadingIndicator = document.getElementById('loadingIndicator');

    e.preventDefault(); // Prevent the default form submission

    var submitBtn = document.getElementById('submitBtn-form');
    currentUrl = window.location.href;
    editing = false
    if (currentUrl.includes('edit')) {
        editing = true
        // pop up that says "are you sure you want to update this reaction?"
        var userConfirmed = confirm('Are you sure you want to update your saved reaction? This action cannot be undone.');
        if (!userConfirmed) {
            return; // Exit the function and do not submit form
        }
    }

    var loadingIndicator = document.getElementById('loadingIndicator');

    var divElement = document.getElementById("organTags");
    var organTags = Array.from(divElement.getElementsByClassName('tag'))
                            .map(tag => tag.firstChild.textContent.trim());

    // Check if the subsystem field is filled
    var subsystemField = document.getElementById('subsystemField').value;
    if (!subsystemField.trim()) {
        alert('Please enter a subsystem.');
        return; // Exit the function and do not submit form
    }

    var inputsGroups = document.querySelectorAll('.inputs-group');
    
    for (var i = 0; i < inputsGroups.length; i++) {
        var group = inputsGroups[i];
        var statusDot = group.querySelector('.status-dot');
        if (statusDot.style.display === 'none') {
            var errorMessage = 'Verify all metabolites before creating reaction.';
            showErrorModal(errorMessage);
            window.scrollTo(0, 0);
            return
        }
    }
    // Continue with form submission

    const disabledInputs = this.querySelectorAll('input:disabled, select:disabled');
    disabledInputs.forEach(input => input.disabled = false);
    var formData = new FormData(this);

    disabledInputs.forEach(input => input.disabled = true);
    var nameData = {};
    var metaboliteFields = document.querySelectorAll('.substrates-name, .products-name');
    var allNamesEntered = true;
    metaboliteFields.forEach(function(input, index) {
        var key = input.name + (index + 1);
        var value = input.value;
        nameData[key] = value;

        if (input.value === '') {
            allNamesEntered = false;
        }
    });

    if (!allNamesEntered) {
        var errorMessage = 'Enter all metabolite names before creating reaction.';
        showErrorModal(errorMessage);
        window.scrollTo(0, 0);
        return; // Exit the function and do not submit form
    }
    submitBtn.disabled = true;
    loadingIndicator.style.display = 'flex';
    if (subsystemList.length === 0) {
        subsystemList = await updateSubsystems();
    }
    var isValidSubsystem = subsystemList.some(subsystem => subsystem.toLowerCase() === subsystemField.toLowerCase());

    if (!isValidSubsystem) {
        var userConfirmed = confirm(`Are you sure you want to add a new subsystem "${subsystemField}"?`);
        if (!userConfirmed) {
            var errorMessage = 'The subsystem entered is not valid.';
            showErrorModal(errorMessage);
            window.scrollTo(0, 0);
            submitBtn.disabled = false;
            loadingIndicator.style.display = 'none';
            return; // Exit the function and do not submit form
        } else {
            if (sessionStorage.getItem('userID') !== null) {
            // Add the new subsystem to the list
            subsystemList.push(subsystemField);
            persistSubsystemList(subsystemField); // Persist the new subsystem
            }
            else {
                var errorMessage = 'Please login to add a new subsystem.';
                showErrorModal(errorMessage);
                window.scrollTo(0, 0);
                submitBtn.disabled = false;
                loadingIndicator.style.display = 'none';
                return; // Exit the function and do not submit form
            }
        }
    }
    var skipAtomMapping = document.getElementById('skipAtomMapping').checked;
    formData.append('skipAtomMapping', skipAtomMapping);
    formData.append('nameData', JSON.stringify(nameData));
    formData.append('organs', JSON.stringify(organTags));
    const url = window.location.href;
    const parsedUrl = new URL(url);
    
    // Extract query parameters
    const params = new URLSearchParams(parsedUrl.search);
    const reactionId = params.get('reaction_id');
    // Get the value of reaction_id
    const action = params.get('action');
    formData.append('action', action);

    if (action === 'edit') {
        formData.append('action', action);

        formData.append('reaction_id', reactionId);
    }   
    formData.append('userID', sessionStorage.getItem('userID'));


    // Check if identical reaction exists
    if (!editing){
        const shouldProceed = await checkIdenticalReaction(formData, loadingIndicator, submitBtn);
        if (!shouldProceed) {
            submitBtn.disabled = false;
            loadingIndicator.style.display = 'none';
            return;
        }
    }
    // Create a new reaction
    fetch(inputReactionUrl, {
        method: 'POST',
        body: formData,
        headers: {
            'X-Requested-With': 'XMLHttpRequest',
            'X-CSRFToken': csrfToken
        }
    })
    .then(response => response.json())
    .then(async (data) => {
        if (data.status === 'error') {
            showErrorModal(data.message);
            window.scrollTo(0, 0);
            submitBtn.disabled = false;
            loadingIndicator.style.display = 'none';
            return;
        }
        else if (action !== 'edit' && data.reaction_id) { // reactionId is taken from the response
            const userID = sessionStorage.getItem('userID');
            const reactionId = data.reaction_id; // Ensure reactionId is obtained from the response
            
            // Send additional request to save CreatedReaction (to keep track of which user created which reaction)
            fetch('create-reaction/', {
                method: 'POST',
                body: JSON.stringify({
                    user_id: userID,
                    reaction_id: reactionId
                }),
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                }
            })
            .then(response => response.json())
            .then(result => {
                if (result.success) {
                    console.log('CreatedReaction saved successfully.');
                } else {
                    console.error('Failed to save CreatedReaction:', result.error);
                }
            })
            .catch(error => {
                console.error('Error saving CreatedReaction:', error);
            });
        }
        
        setTimeout(function() {
            const redirectUrl = window.location.origin + "/?reaction_id=" + data.reaction_id;
            window.location.href = action === 'edit' ? redirectUrl + "&action=edit" : redirectUrl;
        }, 10); 
    })
    .catch(error => {
        console.error('Error:', error);
        const errorMessage = 'An unexpected error occurred.';
        showErrorModal(errorMessage);
        window.scrollTo(0, 0);
        submitBtn.disabled = false;
        loadingIndicator.style.display = 'none';
    });
});    


function persistSubsystemList(newSubsystem) {
    fetch('/update_subsystems/', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrfToken
        },
        body: JSON.stringify({ subsystems: [newSubsystem] })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            console.error('Error in updating subsystems:', data.message);
        } else {
            console.log('Subsystem list updated successfully.');
        }
    })
    .catch(error => console.error('Error:', error));
}


function showLoader() {
    document.getElementById('loadingIndicator').style.display = 'flex';
}

// Function to hide the loading indicator
function hideLoader() {
    document.getElementById('loadingIndicator').style.display = 'none';
}





function showErrorModal(message) {
    var errorMessageElement = document.getElementById('error-message');
    
    // Debugging output

    errorMessageElement.innerText = message;
    document.getElementById('error-message').style.display = 'block';
    document.getElementById('error-modal').style.display = 'block';
}




function updateStatusDots(containerId, foundList, miriamsList) {
    const container = document.getElementById(containerId);
    const statusDots = container.querySelectorAll('.status-dot');
    statusDots.forEach((dot, index) => {
        dot.style.display = 'block'; // Ensure the dot is displayed

        if (foundList[index]) {
            dot.className = 'status-dot found';
            dot.setAttribute('data-tooltip', 'Metabolite found in VMH');
            dot.style.backgroundColor = ''; // Reset color to default
        } else {
            dot.className = 'status-dot not-found';
            dot.setAttribute('data-tooltip', 'Metabolite not found in VMH');
            dot.style.backgroundColor = ''; // Reset color to default
        }

        if (foundList[index] && miriamsList[index]) {
            dot.onclick = () => window.open(miriamsList[index], '_blank');
            dot.style.cursor = 'pointer';
        } else {
            dot.onclick = null;
            dot.style.cursor = 'default';
        }
    });
}

function updateStatusDot(dot, found, miriam) {
    dot.style.display = 'block'; // Ensure the dot is displayed

    if (found) {
        dot.className = 'status-dot found';
        dot.setAttribute('data-tooltip', 'Metabolite found in VMH');
        dot.style.backgroundColor = ''; // Reset color to default
    } else {
        dot.className = 'status-dot not-found';
        dot.setAttribute('data-tooltip', 'Metabolite not found in VMH');
        dot.style.backgroundColor = ''; // Reset color to default
    }

    if (found && miriam) {
        dot.onclick = () => window.open(miriam, '_blank');
        dot.style.cursor = 'pointer';
    } else {
        dot.onclick = null;
        dot.style.cursor = 'default';
    }
}

function confirmAll() {
    // Get all the elements with the class 'done-field-btn-all'
    const verifyAllButtons = document.querySelectorAll('.done-field-btn-all');

    // Loop through each button and trigger a click event
    verifyAllButtons.forEach(button => button.click());
}


function setupTooltips() {
    document.querySelectorAll('.info-symbol').forEach(item => {
        item.addEventListener('mouseenter', function () {
            const tooltipContent = this.getAttribute('data-tooltip-content');
            const tooltip = document.createElement('div');
            tooltip.className = 'tooltip';
            tooltip.innerHTML = tooltipContent;
            this.appendChild(tooltip);
        });
        item.addEventListener('mouseleave', function () {
            this.removeChild(this.querySelector('.tooltip'));
        });
    });
}