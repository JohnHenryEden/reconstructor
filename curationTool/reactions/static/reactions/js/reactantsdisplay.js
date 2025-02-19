document.addEventListener('DOMContentLoaded', async function () {
    $('#savedMetabolitesModal').modal();

    $('#userDropdown').dropdown({
        action: 'hide',
        onChange: function(value, text, $selectedItem) {
            if ($selectedItem.attr('id') === 'dropdown-saved-reactions') {
                userId = sessionStorage.getItem('userID');
                if (userId) {
                    window.location.href = '/saved_reactions';
                }
                else{
                    var errorMessage = 'Please login to view saved reactions.';
                    showErrorModal(errorMessage);
                }
            } else if ($selectedItem.attr('id') === 'dropdown-saved-metabolites') {
                loadSavedMetabolites();
                $('#savedMetabolitesModal').modal('show');
            }
        }
    });
    if (sessionStorage.getItem('userID') !== null) {
        username = sessionStorage.getItem('userName');
        userID = sessionStorage.getItem('userID');
        document.getElementById('userDisplay').innerHTML = `<i class="icon user"></i> User: ${username}`;
        document.getElementById('loginButton').textContent = 'Log out';
        setLoggedInStatusBasedOnUrl('');
        fetch(setSessionUser, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ 'userID': userID })
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    console.log('Session user set successfully:', data.message);
                } else {
                    var errorMessageContainer = 'Error in setting session user: ' + data.message;
                    showErrorModal(errorMessageContainer);   
                }            })

    }

    setupTooltips();
    createnewreaction();
    createGeneInfoInput();
    setupdate();
    displayreactioninfo(reactionData = null);
    attachEventListenersToSelects();
    toggleStructure();
    updateAtomChargeCounters();
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    const action = urlParams.get('action');
    if (reactionId || action === 'edit') {
        fetch(getReaction + reactionId)
            .then(response => {
                if (!response.ok) {
                    console.error(`HTTP error! status: ${response.status}`);
                    throw new Error(`HTTP error! status: ${response.status}`);
                }
                return response.json();
            })
            .then(async reactionData => {
                await updateFormFields(reactionData);
                reactionData.substrates = ['test'];
                confirmAll();
                displayDivs(reactionData);
                if (reactionData.short_name) {
                    let reactionDataName = reactionData.short_name;
                    if(reactionData.short_name.length>20){
                        reactionDataName = reactionData.short_name.substring(0, 30) + "...";
                    }
                    let reactionDataDescription = reactionData.description;
                    reactionStatusInfo = {
                        'name': reactionDataName,
                        'description': reactionDataDescription
                      };
                    setLoggedInStatusBasedOnUrl(reactionStatusInfo);
                    DisplayTag(reactionData.Organs);
                }
            })
            .catch(error => {
                console.error('Error fetching reaction data:', error);
                // Optionally, handle the error by displaying a message to the user
            });
    }
    subsystemList = await updateSubsystems();
    setupTooltips();

});
async function updateSubsystems() {
    if (subsystemList.length === 0) {
        try {
            window.scrollTo(0, 0);

            const response = await fetch(getVMHsubsystems, {
                method: 'GET',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken
                }
            });

            const data = await response.json();

            if (data.error) {
                showErrorModal(data.message);
                return []; // Return an empty list if there's an error
            } else {
                hidemodal();
                return data.subsystem_list; // Return the fetched list
            }
        } catch (error) {
            console.error('Error fetching subsystems:', error);
            return []; // Return an empty list in case of failure
        }
    }
    else {
        return subsystemList;
    }
}
function createnewreaction() {
    // Get the input element by its ID
    const resetbutton = document.getElementById('ResetButton');
    // Add an event listener to the input element to handle the click event
    resetbutton.addEventListener('click', function (event) {
        // Prevent the default form submission behavior
        event.preventDefault();

        // Redirect to the homepage
        window.location.href = window.location.origin;
    });
}

function setupdate(){
    const urlParams = new URLSearchParams(window.location.search);
    const action = urlParams.get('action');
    if (action === 'edit') {
        document.getElementById('submitBtn-form').childNodes[2].nodeValue = 'Update Reaction';
    }

}



