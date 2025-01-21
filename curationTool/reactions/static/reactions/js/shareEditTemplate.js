document.addEventListener('DOMContentLoaded', function () {
  const shareTemplatesContainer = document.getElementById('shareTemplatesContainer');
  const manageTemplateBtn = document.getElementById('manageTemplateBtn');
  const shareTemplateBtn = document.getElementById('shareTemplateBtn');
  const editTemplateBtn = document.getElementById('editTemplateBtn');

  // Open Manage Templates Modal
  manageTemplateBtn.addEventListener('click', function () {
    userID = sessionStorage.getItem('userID');
    if (!userID) {
        alert('Please log in to manage templates.');
        return;
    }
    $('#manageTemplatesModal').modal('show');
  });

  // Open Share Templates Modal when Share button is clicked
  shareTemplateBtn.addEventListener('click', function () {
      $('#manageTemplatesModal').modal('hide'); // Close Manage Templates Modal
      $('#shareTemplateModal').modal('show');  // Open Share Templates Modal
  });

  // Placeholder: Open Edit Templates Modal when Edit button is clicked
  editTemplateBtn.addEventListener('click', function () {
      $('#manageTemplatesModal').modal('hide'); // Close Manage Templates Modal
      $('#editTemplateModal').modal('show');   // Open Edit Templates Modal
  });

  const confirmShareBtn = document.getElementById('confirmShare');
  let templateList = [];
  let templatesFetched = { value: false };

  // Fetch templates and populate buttons
  shareTemplateBtn.addEventListener('click', async function () {
      $('#shareTemplateModal').modal('show');
      shareTemplatesContainer.innerHTML = ''; // Clear existing buttons

      // Fetch templates (reuse fetchTemplates if available)
      await window.ReactionUtils.fetchTemplates(templateList, templatesFetched);

      // Create buttons for each template
      templateList.forEach(templateName => {
          const button = document.createElement('div');
          button.textContent = templateName;
          button.className = 'template-button';
          button.addEventListener('click', function () {
              button.classList.toggle('selected'); // Toggle selection
          });
          shareTemplatesContainer.appendChild(button);
      });
  });

  // Handle share action
  confirmShareBtn.addEventListener('click', async function () {
      const selectedTemplates = [...shareTemplatesContainer.getElementsByClassName('selected')].map(
          button => button.textContent
      );
      const shareWithUser = document.getElementById('shareWithUser').value.trim();
      const userID = sessionStorage.getItem('userID');

      if (!userID) {
          alert('Please log in before sharing templates.');
          return;
      }
      if (!selectedTemplates.length) {
          alert('Please select at least one template to share.');
          return;
      }
      if (!shareWithUser) {
          alert('Please enter the username of the user to share with.');
          return;
      }

      try {
          const response = await fetch('/share_template/', {
              method: 'POST',
              headers: {
                  'Content-Type': 'application/json',
                  'X-Requested-With': 'XMLHttpRequest',
                  'X-CSRFToken': csrfToken,
              },
              body: JSON.stringify({
                  userID: userID,
                  share_with_user: shareWithUser,
                  template_names: selectedTemplates,
              }),
          });

          const data = await response.json();
          if (response.ok && data.status === 'success') {
              alert('Templates shared with ' + shareWithUser + ' successfully');
              $('#shareTemplateModal').modal('hide');
          } else {
              alert(`Error sharing templates: ${data.message || 'Unknown error'}`);
          }
      } catch (error) {
          console.error('Error:', error);
          alert('An error occurred while sharing the templates.');
      }
  });
});

document.addEventListener('DOMContentLoaded', function () {
    const editTemplateBtn = document.getElementById('editTemplateBtn');
    const editTemplatesContainer = document.getElementById('editTemplatesContainer');
    const renameTemplateBtn = document.getElementById('renameTemplateBtn');
    const deleteTemplateBtn = document.getElementById('deleteTemplateBtn');
    const newTemplateName = document.getElementById('newTemplateName');
    userID = sessionStorage.getItem('userID');
    let selectedTemplate = null; // Keep track of the selected template
    let templateList = []; // Template names

    // Open Edit Templates Modal and load templates
    editTemplateBtn.addEventListener('click', async function () {
        $('#editTemplateModal').modal('show');
        editTemplatesContainer.innerHTML = ''; // Clear existing buttons

        // Fetch templates (reuse fetchTemplates if available)
        await window.ReactionUtils.fetchTemplates(templateList, { value: false });

        // Create buttons for each template
        templateList.forEach(templateName => {
            const button = document.createElement('div');
            button.textContent = templateName;
            button.className = 'template-button';
            button.addEventListener('click', function () {
                // Highlight the selected template
                document.querySelectorAll('.template-button').forEach(btn => btn.classList.remove('selected'));
                button.classList.add('selected');
                selectedTemplate = templateName; // Update selected template
            });
            editTemplatesContainer.appendChild(button);
        });
    });

    // Rename Template
    renameTemplateBtn.addEventListener('click', async function () {
        if (!selectedTemplate) {
            alert('Please select a template to rename.');
            return;
        }

        const newName = newTemplateName.value.trim();
        if (!newName) {
            alert('Please enter a new name for the template.');
            return;
        }

        try {
            const response = await fetch('/rename_template/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                },
                body: JSON.stringify({
                    old_name: selectedTemplate,
                    new_name: newName,
                    userID: userID,
                }),
            });

            const data = await response.json();
            if (response.ok && data.status === 'success') {
                alert('Template renamed successfully!');
                editTemplateBtn.click(); // Refresh the template list
                newTemplateName.value = ''; // Clear the input field
            } else {
                alert(`Error renaming template: ${data.message || 'Unknown error'}`);
            }
        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while renaming the template.');
        }
    });

    // Delete Template
    deleteTemplateBtn.addEventListener('click', async function () {
        if (!selectedTemplate) {
            alert('Please select a template to delete.');
            return;
        }

        const confirmed = confirm(`Are you sure you want to delete the template "${selectedTemplate}"?`);
        if (!confirmed) return;

        try {
            const response = await fetch('/delete_template/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                },
                body: JSON.stringify({
                    template_name: selectedTemplate,
                    userID: userID,
                }),
            });

            const data = await response.json();
            if (response.ok && data.status === 'success') {
                alert('Template deleted successfully!');
                editTemplateBtn.click(); // Refresh the list of templates
                newTemplateName.value = ''; // Clear the input field
            } else {
                alert(`Error deleting template: ${data.message || 'Unknown error'}`);
            }
        } catch (error) {
            console.error('Error:', error);
            alert('An error occurred while deleting the template.');
        }
    });
});
