document.addEventListener('DOMContentLoaded', function () {
    const modalButton = document.getElementById('getrxnfromtemplateBtn');
    const reactionField = document.getElementById('reactionField'); 
    const reactionDropdown = document.getElementById('reactionDropdown');
    let reactionList = [];
    let templatesFetched = false; // Track if templates have been fetched

    // Fetch the list of templates
    async function fetchTemplates() {
        if (templatesFetched) return; // Fetch only once
        try {
            userID = sessionStorage.getItem('userID');
            console.log(userID);
            // send the user ID if logged in
            const response = await fetch('/list_templates/', {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ 'userID': userID })
            });

            if (response.ok) {
                const data = await response.json();
                reactionList = data.templates;
                templatesFetched = true; // Mark as fetched
            } else {
                console.error('Failed to fetch templates');
            }
        } catch (error) {
            console.error('Error:', error);
        }
    }

    // Populate dropdown with available templates
    function populateDropdown(list) {
        reactionDropdown.innerHTML = ''; // Clear existing dropdown content
        list.forEach(template => {
            const element = document.createElement('div');
            element.textContent = template;
            element.addEventListener('click', function () {
                reactionField.value = this.textContent; // Set selected value in the input field
                reactionDropdown.style.display = 'none'; // Hide dropdown
            });
            reactionDropdown.appendChild(element);
        });
        reactionDropdown.style.display = list.length > 0 ? 'block' : 'none';
    }

    // Show dropdown on input click
    reactionField.addEventListener('click', function () {
        populateDropdown(reactionList); // Show the full list on click
    });

    // Filter dropdown based on user input
    reactionField.addEventListener('keyup', function () {
        const inputVal = this.value.toLowerCase();
        const matches = reactionList.filter(template => template.toLowerCase().includes(inputVal));
        populateDropdown(matches);
    });

    // Hide dropdown when clicking outside
    document.addEventListener('click', function (event) {
        if (!reactionDropdown.contains(event.target) && !reactionField.contains(event.target)) {
            reactionDropdown.style.display = 'none';
        }
    });

    // Show modal and fetch templates
    modalButton.addEventListener('click', async function () {
        $('#getrxnfromtemplateModal').modal('show');
        await fetchTemplates(); // Fetch templates when modal is opened
    });

    // Apply selected template
    document.getElementById('applyTemplate').addEventListener('click', async function () {
        const selectedValue = reactionField.value;
        if (!reactionList.includes(selectedValue)) {
            alert('Please select a valid template from the dropdown.');
            return;
        }

        try {
            const response = await fetch('/get_rxn_template/', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ reaction_type: selectedValue })
            });
            if (response.ok) {
                const data = await response.json();
                await updateFormFields(data);
                $('.ui.modal.getrxntemplate').modal('hide');
            } else {
                console.error('Failed to fetch template');
            }
        } catch (error) {
            console.error('Error:', error);
        }
    });
});
