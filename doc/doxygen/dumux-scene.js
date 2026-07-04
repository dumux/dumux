// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Renders interactive 3D scenes in the documentation from small .vtp data
// files using vtk.js (loaded from a CDN with a local fallback, see header.html).
//
// Any element with class "dumux-scene" is turned into an interactive viewer.
// Supported data attributes:
//   data-scene (either) name of a geometry registered on window.dumuxScenes by
//              a scene .js file loaded with a <script src>. This works from
//              file:// too, as no data is fetched.
//   data-src   (either) URL of a .vtp file (PolyData) fetched at view time.
//              Only works when the docs are served over http(s).
//   data-field (optional) point data array to color by; defaults to the
//              active scalars, or the first point array
//   data-cmap  (optional) vtk.js color map preset name (default "Cool to Warm")
//
// The scene background and legend follow the current (light/dark) page theme
// and update live when the theme is toggled. The scene data files are produced
// by bin/postprocessing/vtk_to_scene.py.

(function () {
  "use strict";

  // Contexts of all initialized scenes, so the theme can be re-applied on toggle.
  var scenes = [];

  // Wait until vtk.js is available (the CDN copy loads synchronously, but the
  // local fallback in header.html is injected asynchronously on CDN failure).
  function whenVtkReady(callback) {
    if (window.vtk) { callback(); return; }
    var waited = 0;
    var timer = setInterval(function () {
      waited += 100;
      if (window.vtk) { clearInterval(timer); callback(); }
      else if (waited >= 15000) {
        clearInterval(timer);
        console.warn("dumux-scene: vtk.js could not be loaded; 3D scenes disabled.");
      }
    }, 100);
  }

  function parseRgb(str) {
    var m = str && str.match(/rgba?\(([^)]+)\)/);
    if (!m) return null;
    var p = m[1].split(",").map(function (x) { return parseFloat(x); });
    if (p.length < 3) return null;
    if (p.length >= 4 && p[3] === 0) return null; // fully transparent
    return [p[0], p[1], p[2]];
  }

  // Colors derived from the current page theme: the actual page background, a
  // contrasting foreground for the legend, and the page's font family.
  function pageColors() {
    var bg = parseRgb(getComputedStyle(document.body).backgroundColor)
          || parseRgb(getComputedStyle(document.documentElement).backgroundColor)
          || [255, 255, 255];
    var luminance = (0.299 * bg[0] + 0.587 * bg[1] + 0.114 * bg[2]) / 255.0;
    var foreground = luminance < 0.5 ? "#e6e6e6" : "#1a1a1a";
    var fontFamily = getComputedStyle(document.body).fontFamily || "sans-serif";
    return {
      background: [bg[0] / 255.0, bg[1] / 255.0, bg[2] / 255.0],
      foreground: foreground,
      fontFamily: fontFamily,
    };
  }

  function applyTheme(ctx) {
    var colors = pageColors();
    ctx.renderer.setBackground(colors.background);
    if (ctx.scalarBar) {
      var style = {
        fontColor: colors.foreground,
        fontStyle: "normal",
        fontFamily: colors.fontFamily,
      };
      ctx.scalarBar.setAxisTextStyle(Object.assign({ fontSize: 16 }, style));
      ctx.scalarBar.setTickTextStyle(Object.assign({ fontSize: 12 }, style));
    }
    ctx.renderWindow.render();
  }

  function base64ToArrayBuffer(b64) {
    var binary = atob(b64);
    var bytes = new Uint8Array(binary.length);
    for (var i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
    return bytes.buffer;
  }

  // Turn a loaded PolyData into an interactive, theme-aware viewer in container.
  function buildScene(container, polydata, presetName, field) {
    var Core = window.vtk.Rendering.Core;
    var pd = polydata.getPointData();
    var array = field ? pd.getArrayByName(field) : null;
    if (!array) array = pd.getScalars();
    if (!array && pd.getNumberOfArrays() > 0) array = pd.getArrayByIndex(0);
    if (array) pd.setActiveScalars(array.getName());

    var mapper = Core.vtkMapper.newInstance();
    mapper.setInputData(polydata);
    mapper.setInterpolateScalarsBeforeMapping(true);

    var lut = null;
    if (array) {
      var range = array.getRange();
      lut = Core.vtkColorTransferFunction.newInstance();
      var preset = Core.vtkColorTransferFunction.vtkColorMaps.getPresetByName(presetName);
      if (preset) lut.applyColorMap(preset);
      lut.setMappingRange(range[0], range[1]);
      lut.updateRange();
      mapper.setLookupTable(lut);
      mapper.setScalarRange(range[0], range[1]);
      mapper.setColorModeToMapScalars();
      mapper.setScalarModeToUsePointData();
    }

    var actor = Core.vtkActor.newInstance();
    actor.setMapper(mapper);

    var fsrw = window.vtk.Rendering.Misc.vtkFullScreenRenderWindow.newInstance({
      rootContainer: container,
      containerStyle: { width: "100%", height: "100%", position: "relative" },
    });
    var renderer = fsrw.getRenderer();
    renderer.addActor(actor);

    var scalarBar = null;
    if (lut) {
      scalarBar = Core.vtkScalarBarActor.newInstance();
      scalarBar.setScalarsToColors(lut);
      scalarBar.setAxisLabel(array ? array.getName() : "");
      if (scalarBar.setDrawNanAnnotation) scalarBar.setDrawNanAnnotation(false);
      renderer.addActor2D ? renderer.addActor2D(scalarBar) : renderer.addActor(scalarBar);
    }

    renderer.resetCamera();
    var camera = renderer.getActiveCamera();
    camera.azimuth(30);
    camera.elevation(-25);
    renderer.resetCameraClippingRange();

    var ctx = { renderer: renderer, scalarBar: scalarBar, renderWindow: fsrw.getRenderWindow() };
    scenes.push(ctx);
    applyTheme(ctx); // sets theme background + legend style and renders
  }

  function hide(container, message) {
    console.warn("dumux-scene: " + message);
    container.style.display = "none";
  }

  function initScene(container) {
    var field = container.getAttribute("data-field") || "";
    var presetName = container.getAttribute("data-cmap") || "Cool to Warm";
    var reader = window.vtk.IO.XML.vtkXMLPolyDataReader.newInstance();

    // Preferred path: geometry embedded via a <script> tag (works from file://).
    var sceneKey = container.getAttribute("data-scene");
    if (sceneKey) {
      var registry = window.dumuxScenes || {};
      if (!registry[sceneKey]) { hide(container, "scene '" + sceneKey + "' not registered"); return; }
      reader.parseAsArrayBuffer(base64ToArrayBuffer(registry[sceneKey]));
      buildScene(container, reader.getOutputData(0), presetName, field);
      return;
    }

    // Fallback: fetch a .vtp URL (requires the docs to be served over http(s)).
    var src = container.getAttribute("data-src");
    if (!src) { hide(container, "no data-scene or data-src given"); return; }
    reader.setUrl(src).then(function () {
      buildScene(container, reader.getOutputData(0), presetName, field);
    }, function (err) {
      hide(container, "could not load " + src + ": " + err);
    });
  }

  function reapplyTheme() {
    for (var i = 0; i < scenes.length; i++) applyTheme(scenes[i]);
  }

  function watchTheme() {
    // doxygen-awesome toggles dark-mode/light-mode classes on <html>.
    try {
      new MutationObserver(reapplyTheme).observe(
        document.documentElement, { attributes: true, attributeFilter: ["class"] }
      );
    } catch (e) { /* MutationObserver unavailable */ }
    // Also react to the OS color-scheme preference changing.
    try {
      window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", reapplyTheme);
    } catch (e) { /* addEventListener on matchMedia unavailable */ }
  }

  function initAll() {
    var containers = document.querySelectorAll(".dumux-scene");
    for (var i = 0; i < containers.length; i++) initScene(containers[i]);
    watchTheme();
  }

  function onReady() { whenVtkReady(initAll); }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", onReady);
  } else {
    onReady();
  }
})();
