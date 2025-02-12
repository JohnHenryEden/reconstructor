function loadAtomMappingDiv(reactionData) {
    let messageContainer = document.getElementById('reactionFoundMessage');
    let contentDiv = document.querySelector('.content-div[name="atommapping-div"]');
    messageContainer.innerHTML = ''; // Clear previous messages

    if (reactionData.error) {
        let formattedErrorMessage = reactionData.error.replace(/\n/g, '<br>');
        let errorMessage = formattedErrorMessage;
        $('#error-message').html(errorMessage);
        $('#error-modal').modal('show');
        window.scrollTo(0, 0);
    } else {
        let reactionImage = contentDiv.getElementsByClassName('reaction-image')[0];
        let imageClass = reactionImage.className;
        let reactionImageName = reactionImage.name;
        let reactionImageId = reactionImage.id;
        
        // Remove the element and ensure it's done before continuing
        reactionImage.remove();
        
        // Create a new image element with the previous attributes
        let newReactionImage = document.createElement('img');
        newReactionImage.className = imageClass;
        newReactionImage.id = reactionImageId;
        newReactionImage.name = reactionImageName; 
    
        // Generate cache buster and set the new image source
        let cacheBuster = new Date().getTime() + "_" + Math.random();
        newReactionImage.src = MEDIA_URL + reactionData.visualization[0] + '?v=' + cacheBuster; // Assuming data.visualization[0] contains the image path
    
        // Append the new image to the contentDiv
        contentDiv.appendChild(newReactionImage);
    
        let lastMousePos = { x: 0, y: 0 };

        newReactionImage.addEventListener('load', function() {
            let currentScale = 1.5;
        
            // 1) Initialize blowup once.
            $("#" + reactionImageId).blowup({
                width: 300,
                height: 300,
                border: "6px solid #f2711c",
                scale: currentScale
            });
        
            // 2) Track mouse position as we move around
            $("#" + reactionImageId).on("mousemove.blowupPos", function(e) {
                lastMousePos.x = e.pageX;
                lastMousePos.y = e.pageY;
            });
        
            // 3) On wheel scroll, zoom in/out
            $("#" + reactionImageId).on("wheel", function(e) {
                e.preventDefault();
                if (e.originalEvent.deltaY < 0) {
                    currentScale += 0.2; // Zoom in
                } else {
                    currentScale = Math.max(0.2, currentScale - 0.2); // Zoom out
                }
                
                // Remove the old blowup events and lens
                $(this).off("mouseenter mouseleave mousemove");
                $(".BlowupLens").remove();
        
                // Re-initialize blowup with the new scale
                $(this).blowup({
                    width: 300,
                    height: 300,
                    border: "6px solid #f2711c",
                    scale: currentScale
                });
        
                // Re-enable mouse tracking on the new instance
                $(this).on("mousemove.blowupPos", function(e) {
                    lastMousePos.x = e.pageX;
                    lastMousePos.y = e.pageY;
                });
        
                // Force the lens to appear if still hovering
                $(this).trigger("mouseenter");
        
                // 4) Immediately 'fake' a mousemove event at the last known position
                $(this).trigger(
                    $.Event("mousemove", {
                        pageX: lastMousePos.x,
                        pageY: lastMousePos.y
                    })
                );
            });
        });       
    }
}
