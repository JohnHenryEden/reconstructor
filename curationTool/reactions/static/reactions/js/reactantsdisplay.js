document.addEventListener('DOMContentLoaded', async function () {
    
    if (sessionStorage.getItem('userID') !== null) {
        setLoggedInStatusBasedOnUrl();
        username = sessionStorage.getItem('userName');
        userID = sessionStorage.getItem('userID');
        document.getElementById('userDisplay').innerHTML = `<i class="icon user"></i> User: ${username}`;
        document.getElementById('loginButton').textContent = 'Log out';

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
                confirmAll();
                displayDivs(reactionData);
                if (reactionData.short_name) {
                    if(reactionData.short_name.length>20){
                        reactionData = reactionData.short_name.substring(0, 30) + "...";
                        setLoggedInStatusBasedOnUrl(reactionData);
                    }
                    else{setLoggedInStatusBasedOnUrl(reactionData.short_name);}
                }
                DisplayTag(reactionData.Organs);
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



