// get every image, expandButtonToggle and the expand factor
const images = document.querySelectorAll('img');
const expandButton = document.getElementById("expand")
const expand_factor = document.getElementById("expand_scale")
let expansion_scale = 1
let expanded = false;

// create a function to exxpand images when clicked on
function handleImageClick(event){
    const image = event.target
    expansion_scale = expand_factor.value;
    console.log(expansion_scale)
    if(!expanded){
        image.style.transform = `scale(${expansion_scale})`
        image.style.zIndex = 1
        const leftDistance = image.getBoundingClientRect().left;
        if(leftDistance < 0){
            image.style.transform = `translateX(${Math.abs(leftDistance)+10}px) scale(2)`;
        }
        expanded=true
    } else {
        image.style.transform = "scale(1)"
        image.style.position = 'relative'
        image.style.zIndex=0
        expanded=false
    }
}

expandButton.addEventListener('change', ()=>{
    if(expandButton.checked){
        images.forEach(image => {
            image.addEventListener('click', handleImageClick)
        })
        expand_factor.style.display='inline';
    } else {
        images.forEach(image => {
            image.removeEventListener('click', handleImageClick)
            image.style.transform = 'scale(1)'
        })
        expand_factor.style.display = 'none'
    }
})
