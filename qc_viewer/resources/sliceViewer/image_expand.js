// get every image, expandButtonToggle and the expand factor
const images = document.querySelectorAll('img');
const expandButton = document.getElementById("expand")
const expand_factor = document.getElementById("expand_scale")
let expansion_scale = 1
let expanded = false;

// create a function to expand images when clicked on
function handleImageClick(event){
    // get target image and the scaling factor
    const image = event.target
    expansion_scale = expand_factor.value;

    // expand image by the scaling factor if not already expanded
    if(!expanded){
        image.style.transform = `scale(${expansion_scale})`
        image.style.zIndex = 1

        // if it is going to expand off screen shift to the right
        const leftDistance = image.getBoundingClientRect().left;
        if(leftDistance < 0){
            image.style.transform = `translateX(${Math.abs(leftDistance)+10}px) scale(${expansion_scale})`;
        }

        expanded=true

    } else {
        // scale images back to original size
        image.style.transform = "scale(1)"
        image.style.position = 'relative'
        image.style.zIndex=0
        expanded=false
    }
}

// Enables or disables the ability to expand images on click
expandButton.addEventListener('change', ()=>{

    // add listener to enable expansion
    if(expandButton.checked){
        images.forEach(image => {
            image.addEventListener('click', handleImageClick)
        })

        // ensure images expand properly
        expand_factor.style.display='inline';

    } else {

        //remove listener to disable expansion
        images.forEach(image => {
            image.removeEventListener('click', handleImageClick)
            image.style.transform = 'scale(1)'
        })

        expand_factor.style.display = 'none'
    }
})
