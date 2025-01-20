function loadMetaboliteInfoDiv(reactionData) {
    const metaboliteInfoDiv = document.querySelector('.content-div[name="metaboliteinfo-div"]');
    if (metaboliteInfoDiv) {
        metaboliteInfoDiv.innerHTML = `
        <div class="div-header">Metabolite Information</div>
        `;
        fillMetaboliteInfoTab(reactionData);
    } else {
        console.error('Metabolite info div not found.');
    }
}
function fillMetaboliteInfoTab(data) {
    const metaboliteInfoContainer = document.getElementById('metaboliteinfo-div');

    data.metabolite_names.forEach((name, index) => {
        const metaboliteDiv = document.createElement('div');
        metaboliteDiv.classList.add('metabolite');

        const toggleDiv = document.createElement('div');
        toggleDiv.classList.add('metabolite-header');

        const nameElement = document.createElement('h3');
        nameElement.textContent = name;

        const toggleButton = document.createElement('button');
        toggleButton.textContent = 'Show 3D Structure';
        toggleButton.classList.add('toggle-button');
        toggleButton.onclick = function() {
            const structureContainer = this.parentNode.parentNode.querySelector('.structure-container');
            if (structureContainer.style.display === 'none') {
                structureContainer.style.display = 'block';
                this.textContent = 'Hide 3D Structure';
        
                if (!structureContainer.hasAttribute('data-viewer-initialized')) {
                    // Basic setup
                    structureContainer.style.height = '400px';
                    structureContainer.style.width = '400px';
                    structureContainer.style.position = 'relative';
        
                    let config = { backgroundColor: 'white' };
                    let viewer = new $3Dmol.createViewer(structureContainer, config);
                    let molecularData = data.metabolite_mol_file_strings[index];
                    let format = 'sdf';
        
                    // Add the molecule model
                    let model = viewer.addModel(molecularData, format);
        
                    viewer.setStyle({}, { stick: { colorscheme: 'greenCarbon' } });
                    let stereoLocations = data.stereo_locations_list[index];
                    for (const loc of stereoLocations) {
                        viewer.setStyle({ model: model, index: loc }, 
                                        { stick: {color: 'magenta' } });
                    }
        
                    // Make atoms clickable if you want labels on click
                    viewer.setClickable({}, true, function(atom, _viewer, _event, _container) {
                        viewer.addLabel(atom.atom, { 
                            position: atom, 
                            backgroundColor: 'darkgreen', 
                            backgroundOpacity: 0.8 
                        });
                    });
        
                    // Zoom and render
                    viewer.zoomTo();
                    viewer.render();
        
                    // Mark as initialized so we don't rebuild next time
                    structureContainer.setAttribute('data-viewer-initialized', 'true');
                }
            } else {
                structureContainer.style.display = 'none';
                this.textContent = 'Show 3D Structure';
            }
        };
        
        toggleDiv.appendChild(nameElement);
        toggleDiv.appendChild(toggleButton);
        metaboliteDiv.appendChild(toggleDiv);
        
        // Existing formula display
        const formulaElement = document.createElement('p');
        formulaElement.textContent = `Charged Formula: ${data.metabolite_formulas[index]}`;
        metaboliteDiv.appendChild(formulaElement);
        
        // NEW: Stereo info simplified
        const stereoCount = data.stereo_counts[index];
        const stereoCountElement = document.createElement('p');
        if (stereoCount > 0) {
            stereoCountElement.textContent = `Number of unspecified Stereo Centers: ${stereoCount} (magenta in the 3D viewer)`;
        } else {
            stereoCountElement.textContent = 'No unspecified stereo centers detected.';
        }
        metaboliteDiv.appendChild(stereoCountElement);
        
        // Structure container
        const structureContainer = document.createElement('div');
        structureContainer.className = 'structure-container';
        structureContainer.style.display = 'none';
        metaboliteDiv.appendChild(structureContainer);
        
        metaboliteInfoContainer.appendChild(metaboliteDiv);
        
    });
}


function toggleStructure() {
    const buttons = document.querySelectorAll('.toggle-button');
    buttons.forEach(button => {
        button.addEventListener('click', () => {
            const messageContainer = button.parentNode.nextElementSibling;
            if (messageContainer.style.display === 'none' || messageContainer.style.display === '') {
                messageContainer.style.display = 'block';
                button.textContent = 'Hide Message';
            } else {
                messageContainer.style.display = 'none';
                button.textContent = 'Show Message';
            }
        });
    });
}
