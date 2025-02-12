function loadChemInfoDiv(reactionData) {
    const contentDiv = document.querySelector('div.content-div[name="cheminfo-div"]');
    contentDiv.innerHTML = `<div class="div-header">Chemical Information</div>`;
    const maindiv = document.createElement('div');

    // Add balanced text
    const balancedText = createBalancedText(reactionData);
    maindiv.appendChild(balancedText);

    // Add molecular formula
    const molecularFormula = createMolecularFormula(reactionData);
    maindiv.appendChild(molecularFormula);

    // Add mass balance information
    const balancedInfoDiv = createInfoDiv(
        `Mass Balanced? ${reactionData.balanced_count[0] ? 'Yes' : 'No'}`,
        'balanced-info'
    );
    maindiv.appendChild(balancedInfoDiv);

    // Add charge balance information
    const chargeBalancedInfoDiv = createInfoDiv(
        `Charge Balanced? ${reactionData.balanced_charge[0] ? 'Yes' : 'No'}`,
        'balanced-info'
    );
    maindiv.appendChild(chargeBalancedInfoDiv);

    // Add atom comparison table
    const comparisonListDiv = document.createElement('div');
    comparisonListDiv.className = 'atoms-comparison';
    populateAtomComparisonTable(
        reactionData.subs_atoms[0],
        reactionData.prods_atoms[0],
        comparisonListDiv,
        { substrates: reactionData.subs_charge[0], products: reactionData.prods_charge[0] },
        reactionData.symb_to_name[0]
    );
    maindiv.appendChild(comparisonListDiv);



    maindiv.style.textAlign = 'center';
    contentDiv.appendChild(maindiv);

    // Add VMH formula section
    appendVMHFormulaSection(contentDiv, reactionData);
}

// Helper function to create the balanced text
function createBalancedText(reactionData) {
    const balancedText = document.createElement('p');
    const isBalanced = reactionData.balanced_count[0] && reactionData.balanced_charge[0];
    balancedText.textContent = isBalanced ? 'Balanced' : 'Not Balanced';
    balancedText.className = `balanced-text ${isBalanced ? 'balanced' : 'not-balanced'}`;
    return balancedText;
}

// Helper function to create molecular formula
function createMolecularFormula(reactionData) {
    const formulaParagraph = document.createElement('p');
    formulaParagraph.textContent = 'Molecular Formula: ';
    const atomColorMap = generateAtomColorMap(reactionData);
    const formulaSpan = document.createElement('span');
    formulaSpan.innerHTML = reactionData.molc_formula[0].replace(/([A-Z][a-z]?)/g, (match) => {
        return atomColorMap[match] ? `<span style="color: ${atomColorMap[match]}">${match}</span>` : match;
    });
    formulaParagraph.appendChild(formulaSpan);
    return formulaParagraph;
}

// Helper function to create generic information divs
function createInfoDiv(text, className) {
    const infoDiv = document.createElement('div');
    infoDiv.className = className;
    infoDiv.textContent = text;
    return infoDiv;
}

// Helper function to generate atom color map
function generateAtomColorMap(reactionData) {
    const atomTypes = new Set([
        ...Object.keys(reactionData.subs_atoms[0]),
        ...Object.keys(reactionData.prods_atoms[0])
    ]);
    const colors = getDistinctColors(atomTypes.size);
    const atomColorMap = {};
    let i = 0;
    atomTypes.forEach((atom) => {
        atomColorMap[atom] = colors[i++];
    });
    return atomColorMap;
}


function getDistinctColors(count) {
    // Function to generate distinct colors for atom types
    const colors = [];
    for (let i = 0; i < count; i++) {
        const hue = (i * 137.508) % 360; // Use golden angle approximation
        colors.push(`hsl(${hue}, 100%, 50%)`);
    }
    return colors;
}

function appendVMHFormulaSection(contentDiv, reactionData) {
    const formulaSection = document.createElement('div');
    formulaSection.className = 'vmh-formula-section';
    formulaSection.style.marginTop = '30px';
    formulaSection.style.padding = '20px';
    formulaSection.style.border = '1px solid #ccc';
    formulaSection.style.borderRadius = '8px';
    formulaSection.style.backgroundColor = '#f4f4f4';
    formulaSection.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.1)';
    formulaSection.style.textAlign = 'center';

    const formulaTitle = document.createElement('h3');
    formulaTitle.textContent = 'VMH Formula';
    formulaTitle.style.marginBottom = '20px';
    formulaSection.appendChild(formulaTitle);

    const formulaDisplay = document.createElement('div');
    formulaDisplay.className = 'vmh-formula-display';
    formulaDisplay.style.fontSize = '18px';
    formulaDisplay.style.fontFamily = 'monospace';
    formulaDisplay.style.padding = '10px';
    formulaDisplay.style.border = '1px dashed #ccc';
    formulaDisplay.style.borderRadius = '4px';
    formulaDisplay.style.backgroundColor = '#fff';
    formulaDisplay.style.cursor = 'pointer';

    if (reactionData.rxn_formula) {
        formulaDisplay.textContent = reactionData.rxn_formula;
    } else {
        const formulaText = [];
        reactionData.subs_sch.forEach((stoich, index) => {
            if (stoich > 1) formulaText.push(parseFloat(stoich).toFixed(1) + ' ');
            formulaText.push(
                reactionData.subs_found[index]
                    ? `${reactionData.substrates[index]}[${reactionData.subs_comps[index]}]` 
                    : `${reactionData.substrates_names[index]}[${reactionData.subs_comps[index]}]`
            );
            if (index < reactionData.subs_sch.length - 1) formulaText.push(' + ');
        });
        
        formulaText.push(
            reactionData.direction === 'forward'
                ? ' -> '
                : reactionData.direction === 'bidirectional'
                ? ' <=> '
                : ''
        );
        
        reactionData.prod_sch.forEach((stoich, index) => {
            if (stoich > 1) formulaText.push(parseFloat(stoich).toFixed(1) + ' ');
            formulaText.push(
                reactionData.prod_found[index]
                    ? `${reactionData.products[index]}[${reactionData.prods_comps[index]}]` 
                    : `${reactionData.products_names[index]}[${reactionData.prods_comps[index]}]`
            );
            if (index < reactionData.prod_sch.length - 1) formulaText.push(' + ');
        });
        
        formulaDisplay.textContent = formulaText.join('');
        
    }

    formulaDisplay.addEventListener('click', () => {
        const text = formulaDisplay.textContent;
        if (navigator.clipboard && navigator.clipboard.writeText) {
            // Use the Clipboard API if available
            navigator.clipboard.writeText(text).then(() => {
                showToast('Formula copied to clipboard!');
            }).catch(err => {
                console.error('Could not copy text: ', err);
            });
        } else {
            // Fallback for browsers that do not support the Clipboard API
            const textarea = document.createElement('textarea');
            textarea.value = text;
            document.body.appendChild(textarea);
            textarea.select();
            try {
                document.execCommand('copy');
                showToast('Formula copied to clipboard!');
            } catch (err) {
                console.error('Could not copy text: ', err);
            }
            document.body.removeChild(textarea);
        }
    });
        
    formulaSection.appendChild(formulaDisplay);

    const generateAbbrButton = document.createElement('button');
    generateAbbrButton.textContent = 'Generate Abbreviations';
    generateAbbrButton.style.marginTop = '20px';
    generateAbbrButton.style.padding = '10px 20px';
    generateAbbrButton.style.fontSize = '16px';
    generateAbbrButton.style.borderRadius = '4px';
    generateAbbrButton.style.border = 'none';
    generateAbbrButton.style.backgroundColor = '#007BFF';
    generateAbbrButton.style.color = '#fff';
    generateAbbrButton.style.cursor = 'pointer';

    generateAbbrButton.addEventListener('click', () => {
        openAbbreviationModal(reactionData);
    });

    formulaSection.appendChild(generateAbbrButton);
    contentDiv.appendChild(formulaSection);
}

function openAbbreviationModal(reactionData) {
    // Create the modal container
    const modal = document.createElement('div');
    modal.style.position = 'fixed';
    modal.style.top = '50%';
    modal.style.left = '50%';
    modal.style.transform = 'translate(-50%, -50%)';
    modal.style.zIndex = '1000';
    modal.style.backgroundColor = '#fff';
    modal.style.padding = '20px';
    modal.style.borderRadius = '8px';
    modal.style.boxShadow = '0 4px 8px rgba(0, 0, 0, 0.2)';
    modal.style.width = '400px';
    modal.style.textAlign = 'center';

    // Modal title
    const modalTitle = document.createElement('h4');
    modalTitle.textContent = 'Generate Abbreviation for Non-VMH Metabolites';
    modalTitle.style.marginBottom = '20px';
    modal.appendChild(modalTitle);

    // Create the loader element
    const loader = document.createElement('div');
    loader.style.display = 'none';
    loader.style.position = 'absolute';
    loader.style.top = '50%';
    loader.style.left = '50%';
    loader.style.transform = 'translate(-50%, -50%)';
    loader.style.padding = '20px';
    loader.style.backgroundColor = 'white';
    loader.style.color = 'black';
    loader.style.borderRadius = '5px';
    loader.style.zIndex = '1000';
    loader.style.textAlign = 'center';
    loader.style.boxShadow = '0px 10px 15px rgba(0, 0, 0, 0.3)';
    loader.textContent = 'Loading...';
    modal.appendChild(loader);
// Collect unique metabolites not found in VMH
const uniqueNotFoundMetabolites = new Map();

reactionData.substrates_names.forEach((name, index) => {
    if (!reactionData.subs_found[index]) {
        uniqueNotFoundMetabolites.set(name, {
            metabolite: reactionData.substrates[index],
            mtype: reactionData.subs_types[index],
        });
    }
});

reactionData.products_names.forEach((name, index) => {
    if (!reactionData.prod_found[index]) {
        // Avoid overwriting if the metabolite is already in the map
        if (!uniqueNotFoundMetabolites.has(name)) {
            uniqueNotFoundMetabolites.set(name, {
                metabolite: reactionData.products[index],
                mtype: reactionData.prods_types[index],
            });
        }
    }
});
// Add close button
const closeButton = document.createElement('button');
closeButton.textContent = 'Close';
closeButton.style.marginTop = '20px';
closeButton.style.padding = '10px 20px';
closeButton.style.fontSize = '16px';
closeButton.style.borderRadius = '4px';
closeButton.style.border = 'none';
closeButton.style.backgroundColor = '#dc3545';
closeButton.style.color = '#fff';
closeButton.style.cursor = 'pointer';

closeButton.addEventListener('click', () => {
    modal.remove();
});

// Generate abbreviations for unique metabolites
uniqueNotFoundMetabolites.forEach((details, name) => {
    const metaboliteDiv = document.createElement('div');
    metaboliteDiv.style.marginBottom = '10px';

    const nameSpan = document.createElement('span');
    nameSpan.textContent = name;
    metaboliteDiv.appendChild(nameSpan);

    const generateButton = document.createElement('button');
    generateButton.textContent = 'Generate';
    generateButton.style.marginLeft = '10px';
    generateButton.style.padding = '5px 10px';
    generateButton.style.fontSize = '14px';
    generateButton.style.borderRadius = '4px';
    generateButton.style.border = 'none';
    generateButton.style.backgroundColor = '#28a745';
    generateButton.style.color = '#fff';
    generateButton.style.cursor = 'pointer';
    function escapeRegExp(string) {
        return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      }      

    generateButton.addEventListener('click', () => {
        showLoaderdiv(generateButton, closeButton); // Pass the closeButton
        fetch('/create-formula-abbr/', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                metabolite: details.metabolite,
                mtype: details.mtype,
                metabolite_name: name,
            }),
        })
            .then((response) => response.json())
            .then((data) => {
                console.log(data);
                if (data.abbr) {
                    hideLoaderdiv(generateButton, closeButton); // Pass the closeButton
    
                    // Update the name span to show the abbreviation
                    nameSpan.textContent = data.abbr;
    
                    const formulaDisplay = document.querySelector('.vmh-formula-display');
    
                    if (formulaDisplay) {
                        const escapedName = escapeRegExp(name);
                        const regex = new RegExp(`\\b${escapedName}\\[(c|g|e|l|m|n|x|r)\\]`, 'gi');
    
                        // Replace the matched metabolites with the abbreviation and the same compartment
                        const updatedFormula = formulaDisplay.textContent.replace(
                            regex,
                            `${data.abbr}[$1]`
                        );
    
                        formulaDisplay.textContent = updatedFormula;
                    }
                    // Save the updated formula to the backend
                    const updatedFormulaText = document.querySelector('.vmh-formula-display').textContent;
                    const reactionId = reactionData.reaction_id; // Assume the reaction ID is available
                    fetch('/save-formula/', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json',
                            'X-CSRFToken': csrfToken, // Replace csrfToken with the actual CSRF token variable
                        },
                        body: JSON.stringify({
                            formula: updatedFormulaText,
                            reaction_id: reactionId,
                        }),
                    })
                        .then((saveResponse) => {
                            if (saveResponse.ok) {
                                console.log('Formula successfully saved.');
                            } else {
                                console.error('Failed to save the formula.');
                            }
                        })
                        .catch((saveError) => console.error('Error saving formula:', saveError));
    
                    // Remove the generate button
                    generateButton.remove();
                }
            })
            .catch((error) => {
                console.error('Error:', error);
                hideLoaderdiv(generateButton, closeButton); 
            });
    });
    
    metaboliteDiv.appendChild(generateButton);
    modal.appendChild(metaboliteDiv);
    });

    modal.appendChild(closeButton);

    // Append modal to document body
    document.body.appendChild(modal);
    
}


function populateAtomComparisonTable(substratesData, productsData, listDiv, charge, symb_to_name) {
    listDiv.style.display = 'flex';
    listDiv.style.justifyContent = 'space-between';
    listDiv.style.alignItems = 'flex-start';
    listDiv.style.padding = '20px';
    listDiv.style.border = '1px solid #ddd';
    listDiv.style.borderRadius = '8px';
    listDiv.style.marginTop = '5%';
    listDiv.style.backgroundColor = '#f9f9f9';
    listDiv.style.fontSize = '16px';  // Replace '14px' with your desired font size

    const allAtoms = { ...substratesData, ...productsData };

    const table = document.createElement('table');
    table.className = 'atom-comparison-table';
    table.style.width = '70%';
    table.style.borderCollapse = 'collapse';

    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');
    ['Atom', 'Substrates', 'Products'].forEach(text => {
        const th = document.createElement('th');
        th.textContent = text;
        th.style.borderBottom = '2px solid #000';
        th.style.padding = '10px';
        th.style.backgroundColor = '#eaeaea';
        th.style.textAlign = 'center';
        th.style.fontWeight = 'bold';
        th.style.fontSize = '19px';  // Replace '16px' with your desired font size
        headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);

    const tbody = document.createElement('tbody');
    Object.keys(allAtoms).forEach(atom => {
        const row = document.createElement('tr');
        const substratesCount = substratesData[atom] || 0;
        const productsCount = productsData[atom] || 0;

        if (substratesCount !== productsCount) {
            row.style.backgroundColor = '#ffcccc'; // Highlight in red
        }

        [symb_to_name[atom] || atom, substratesCount, productsCount].forEach((text, index) => {
            const cell = document.createElement('td');
            cell.textContent = text;
            cell.style.padding = '8px';
            cell.style.borderBottom = '1px solid #ddd';
            cell.style.textAlign = 'center';
            if (index === 0) {
                cell.style.fontWeight = 'bold';
            }
            row.appendChild(cell);
        });
        tbody.appendChild(row);
    });
    table.appendChild(tbody);

    listDiv.appendChild(table);

    const chargeDiv = document.createElement('div');
    chargeDiv.innerHTML = `<strong>Charge Info</strong><br>Substrates: ${charge.substrates}<br>Products: ${charge.products}`;
    chargeDiv.style.marginLeft = '20px';
    chargeDiv.style.padding = '10px';
    chargeDiv.style.border = '1px solid #ccc';
    chargeDiv.style.borderRadius = '8px';
    chargeDiv.style.backgroundColor = '#fff';
    chargeDiv.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
    listDiv.appendChild(chargeDiv);
}
function showLoaderdiv(button, closeButton) {
    button.disabled = true; // Disable the button
    button.style.backgroundColor = '#ccc'; 
    button.style.cursor = 'not-allowed'; 
    button.textContent = 'Generating...'; 
    // Add a spinner and loading message over the button
    const spinner = document.createElement('span');
    spinner.classList.add('spinner');
    spinner.style.marginLeft = '10px';
    spinner.style.display = 'inline-block'; // Ensure it is inline with the button text
    spinner.style.border = '3px solid rgba(0, 0, 0, 0.1)'; // Light gray border
    spinner.style.borderTop = '3px solid #3498db'; // Contrasting blue for spinning effect
    spinner.style.borderRadius = '50%';
    spinner.style.width = '16px';
    spinner.style.height = '16px';
    spinner.style.animation = 'spin 1s linear infinite';
    button.appendChild(spinner);    

    // Disable close button
    closeButton.disabled = true;
    closeButton.style.backgroundColor = '#aaa';
    closeButton.style.cursor = 'not-allowed';
}

function hideLoaderdiv(button, closeButton) {
    button.disabled = false; // Re-enable the button
    button.style.backgroundColor = '#28a745'; // Restore original color
    button.style.cursor = 'pointer'; // Restore cursor

    // Remove the spinner
    const spinner = button.querySelector('.spinner');
    if (spinner) {
        spinner.remove();
    }

    // Re-enable close button
    closeButton.disabled = false;
    closeButton.style.backgroundColor = '#dc3545';
    closeButton.style.cursor = 'pointer';
}
