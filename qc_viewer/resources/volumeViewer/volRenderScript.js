import * as THREE from 'three';

import { GUI } from 'three/addons/libs/lil-gui.module.min.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { VolumeRenderShader1 } from 'three/addons/shaders/VolumeShader.js';
import { VRButton } from 'three/addons/webxr/VRButton.js';

/* 
Script is a reworked version of a threejs example
refined for array data stored in json format 

Copyright © 2010-2024 three.js authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE
*/

let renderer,
    scene,
    camera,
    controls,
    material,
    volconfig,
    cmtextures,
    gui;

let mesh;

init();

function init() {

    scene = new THREE.Scene();
    // Create renderer
    renderer = new THREE.WebGLRenderer();
    renderer.setPixelRatio( window.devicePixelRatio );
    renderer.setSize( window.innerWidth, window.innerHeight );
    renderer.xr.enabled = true;
    document.body.appendChild( renderer.domElement );
    document.body.appendChild(VRButton.createButton(renderer))

    // Create camera (The volume renderer does not work very well with perspective yet)
    const h = 512; // frustum height
    const aspect = window.innerWidth / window.innerHeight;
    camera = new THREE.OrthographicCamera( - h * aspect / 2, h * aspect / 2, h / 2, - h / 2, 1, 1000 );
    camera.position.set( - 64, - 64, 128 );
    camera.up.set( 0, 0, 1 ); // In our data, z is up

    // Create controls
    controls = new OrbitControls( camera, renderer.domElement );
    controls.addEventListener( 'change', render );
    controls.target.set( 64, 64, 128 );			
    controls.minZoom = 0.5;
    controls.maxZoom = 4;
    controls.enablePan = false;
    controls.update();
    


    // The gui for interaction
    volconfig = {channel: 0, cmax:30000, cmin: 500, clim1: 0, clim2: 1, renderstyle: 'iso', isothreshold: 0.05, colormap: 'viridis', channel: 0, zmax: 0,zmin: 0, ymax: 0, ymin: 0, xmax:0, xmin: 0 };
    gui = new GUI();
    gui.add( volconfig, 'clim1', 0, 1, 0.01 ).onChange();
    gui.add( volconfig, 'clim2', 0, 1, 0.01 ).onChange( updateUniforms );
    gui.add( volconfig, 'colormap', { gray: 'gray', viridis: 'viridis' } ).onChange( updateUniforms );
    gui.add( volconfig, 'renderstyle', { mip: 'mip', iso: 'iso' } ).onChange( updateUniforms );
    gui.add( volconfig, 'isothreshold', 0, 1, 0.01 ).onChange( updateUniforms );

    // Load the 4D array from the json file produced
    load(false)

    window.addEventListener( 'resize', onWindowResize );

}

function load(refresh){
    new THREE.FileLoader().load( 'volumeData.json', function ( volume ) {
        volume = JSON.parse(volume)[volconfig.channel]
        // get the length for each axes x,y,z; will only process one channel
        let z_length = volume.length;
        let y_length = volume[0].length;
        let x_length = volume[0][0].length;
        if(!refresh){
            volconfig.zmax = z_length;
            volconfig.ymax = y_length;
            volconfig.xmax = x_length;
            gui.add(volconfig, 'channel', 0, 1, 1).onFinishChange(()=>load(true));
            gui.add(volconfig, 'cmax', 0, 30000, 10).onFinishChange(()=>load(true));
            gui.add(volconfig, 'cmin', 0, 30000, 10).onFinishChange(()=>load(true));
            gui.add(volconfig, 'zmax', 1, z_length, 1 ).onFinishChange( ()=>load(true) );
            gui.add(volconfig, 'zmin', 0, z_length, 1 ).onFinishChange( ()=>load(true) )
            gui.add(volconfig, 'ymax', 1, y_length, 1 ).onFinishChange( ()=>load(true) );
            gui.add(volconfig, 'ymin', 0, y_length, 1 ).onFinishChange( ()=>load(true) );
            gui.add(volconfig, 'xmax', 1, x_length, 1).onFinishChange(()=>load(true));
            gui.add(volconfig, 'xmin', 0, x_length, 1 ).onFinishChange( ()=>load(true) );
        } else {
            scene.remove(mesh)
        }
        // create a new array to transform the array to 1D
        let newData = new Float32Array(x_length*y_length*z_length);

        // loop through every data point in the array
        for(let z=volconfig.zmin; z<volconfig.zmax; z++){
            for(let y = volconfig.ymin; y<volconfig.ymax; y++){
                for(let x = volconfig.xmin; x< volconfig.xmax; x++){
                    // process the point of data an dinput it into the array
                    let zslice = volume.length-1 - z;
                    let dp = volume[zslice][y][x]
                    if(true){
                        let index = x + y*x_length + z*x_length*y_length;
                        dp = (dp-volconfig.cmin)/(volconfig.cmax-volconfig.cmin);
                        newData[index] = dp;
                    }
                }
            }
        }


        // create an object to be processed
        volume = {data: newData, xLength : x_length, yLength : y_length, zLength : z_length};
        // Texture to hold the volume. We have scalars, so we put our data in the red channel.
        // THREEJS will select R32F (33326) based on the THREE.RedFormat and THREE.FloatType.
        const texture = new THREE.Data3DTexture( volume.data, volume.xLength, volume.yLength, volume.zLength );
        texture.format = THREE.RedFormat;
        texture.type = THREE.FloatType;
        texture.minFilter = texture.magFilter = THREE.LinearFilter;
        texture.unpackAlignment = 1;
        texture.needsUpdate = true;

        // Colormap textures
        cmtextures = {
            viridis: new THREE.TextureLoader().load( 'cm_viridis.png', render ),
            gray: new THREE.TextureLoader().load( 'cm_gray.png', render )
        };

        // Material
        const shader = VolumeRenderShader1;

        const uniforms = THREE.UniformsUtils.clone( shader.uniforms );

        uniforms[ 'u_data' ].value = texture;
        uniforms[ 'u_size' ].value.set( volume.xLength, volume.yLength, volume.zLength );
        uniforms[ 'u_clim' ].value.set( volconfig.clim1, volconfig.clim2 );
        uniforms[ 'u_renderstyle' ].value = volconfig.renderstyle == 'mip' ? 0 : 1; // 0: MIP, 1: ISO
        uniforms[ 'u_renderthreshold' ].value = volconfig.isothreshold; // For ISO renderstyle
        uniforms[ 'u_cmdata' ].value = cmtextures[ volconfig.colormap ];

        material = new THREE.ShaderMaterial( {
            uniforms: uniforms,
            vertexShader: shader.vertexShader,
            fragmentShader: shader.fragmentShader,
            side: THREE.BackSide // The volume shader uses the backface as its "reference point"
        } );

        // THREE.Mesh
        const geometry = new THREE.BoxGeometry( volume.xLength, volume.yLength, volume.zLength);
        geometry.translate( volume.xLength / 2 - 0.5, volume.yLength / 2 -0.5, volume.zLength / 2 - 0.5 );

        mesh = new THREE.Mesh( geometry, material );
        scene.add( mesh );
        refresh=true;
        render();
    } );

}

function updateUniforms() {

    material.uniforms[ 'u_clim' ].value.set( volconfig.clim1, volconfig.clim2 );
    material.uniforms[ 'u_renderstyle' ].value = volconfig.renderstyle == 'mip' ? 0 : 1; // 0: MIP, 1: ISO
    material.uniforms[ 'u_renderthreshold' ].value = volconfig.isothreshold; // For ISO renderstyle
    material.uniforms[ 'u_cmdata' ].value = cmtextures[ volconfig.colormap ];

    render();

}

function onWindowResize() {

    renderer.setSize( window.innerWidth, window.innerHeight );

    const aspect = window.innerWidth / window.innerHeight;

    const frustumHeight = camera.top - camera.bottom;

    camera.left = - frustumHeight * aspect / 2;
    camera.right = frustumHeight * aspect / 2;

    camera.updateProjectionMatrix();

    render();

}

function render() {

    renderer.render( scene, camera );
}
