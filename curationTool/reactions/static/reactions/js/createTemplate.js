document.getElementById('createTemplateBtn').addEventListener('click', function() {
    // check if the user is logged in
    if (!sessionStorage.getItem('userID')) {
        alert('Please log in to create a template.');
        return;
    }
    // Hide the first modal
    $('#getrxnfromtemplateModal').modal('hide');

    // Enter "Create Template" mode
    toggleCreateTemplateMode();
});

function toggleCreateTemplateMode() {
    const statusElement = document.getElementById('statusTitle');
    const dotClass = 'dot-red';
    const saveButton = document.getElementById('saveReactionButton');
    let createTemplateButton = document.getElementById('createTemplateButton');
    const resetButton = document.getElementById('ResetButton');
    const formButtonsContainer = document.getElementById('FormButtons');

    // Update status
    statusElement.innerHTML = `Status: <span class="status-dot-top ${dotClass}"></span> Creating Template`;

    // Hide save button
    saveButton.style.display = 'none';
    // Hide create reaction button
    document.getElementById('submitBtn-form').style.display = 'none';
    // Hide skip-atom-mapping
    document.getElementsByClassName('skip-atom-mapping')[0].style.display = 'none';
    // Hide getbuttoncontainer class
    document.getElementsByClassName('getbuttoncontainer')[0].style.display = 'none';

    // Create and show 'Create Template' button if it doesn't exist
    if (!createTemplateButton) {
        createTemplateButton = document.createElement('button');
        createTemplateButton.id = 'createTemplateButton';
        createTemplateButton.className = 'ui button';
        createTemplateButton.textContent = 'Create Template';

        // Copy styles from submitBtn-form
        const submitButton = document.getElementById('submitBtn-form');
        if (submitButton) {
            createTemplateButton.style.color = submitButton.style.color;
            createTemplateButton.innerHTML = submitButton.innerHTML.replace(
                'Create Reaction', 
                'Create Template'
            );
        }

        // Place 'Create Template' button before the ResetButton
        formButtonsContainer.insertBefore(createTemplateButton, resetButton);

        // Add event listener for creating the template
        createTemplateButton.addEventListener('click', function() {
            // Collect form data and send it to the backend
            createTemplate();
        });
    } else {
        createTemplateButton.style.display = 'inline-block';
    }
}

async function createTemplate() {
    const reactionForm = document.getElementById('reactionForm');
    const formData = new FormData(reactionForm);

    // Prompt user for the template name
    const templateName = prompt('Enter a name for the template:');
    if (!templateName) {
        alert('Template name is required.');
        return;
    }

    formData.append('template_name', templateName);

    // Collect organ tags
    const organTagsDiv = document.getElementById('organTags');
    const organTags = Array.from(organTagsDiv.getElementsByClassName('tag'))
        .map(tag => tag.firstChild.textContent.trim());
    formData.append('organs', JSON.stringify(organTags));

    // Validate subsystem field
    const subsystemField = document.getElementById('subsystemField').value.trim();
    if (!subsystemField) {
        alert('Please enter a subsystem.');
        return;
    }

    const isValidSubsystem = subsystemList.some(subsystem => subsystem.toLowerCase() === subsystemField.toLowerCase());
    if (!isValidSubsystem) {
        const userConfirmed = confirm(`Are you sure you want to add a new subsystem "${subsystemField}"?`);
        if (!userConfirmed) {
            alert('The subsystem entered is not valid.');
            return;
        } else if (!sessionStorage.getItem('userID')) {
            alert('Please log in to add a new subsystem.');
            return;
        } else {
            subsystemList.push(subsystemField); // Add to local list
            persistSubsystemList(subsystemField); // Save subsystem if logged in
        }
    }
    formData.append('subsystem', subsystemField);

    // Add the user ID if logged in
    const userID = sessionStorage.getItem('userID');
    if (userID) {
        formData.append('userID', userID);
    }

    // Send the template creation request
    try {
        const response = await fetch('/create_template/', {
            method: 'POST',
            body: formData,
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': document.querySelector('[name=csrfmiddlewaretoken]').value
            }
        });

        if (response.ok) {
            const data = await response.json();
            alert('Template created successfully!');
            // refresh the page
            location.reload();
        } else {
            const errorData = await response.json();
            alert(`Error creating template: ${errorData.message}`);
        }
    } catch (error) {
        console.error('Error creating template:', error);
        alert('An error occurred while creating the template.');
    }
}
