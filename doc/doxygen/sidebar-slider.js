// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later

(function () {
    const root = document.documentElement;
    const storageKey = "sideNavWidth";

    const minWidth = 240;
    const maxWidth = 760;

    function getCurrentWidth() {
        return (
            parseInt(
                getComputedStyle(root).getPropertyValue("--side-nav-fixed-width"),
                10
            ) || minWidth
        );
    }

    function applyWidth(px) {
        const width = Math.max(minWidth, Math.min(maxWidth, px));
        root.style.setProperty("--side-nav-fixed-width", width + "px");
        localStorage.setItem(storageKey, String(width));
    }

    function ensureResizeHandle() {
        const sideNav = document.getElementById("side-nav");
        if (!sideNav) return;

        let handle = document.getElementById("sidebar-resize-handle");

        if (!handle) {
            handle = document.createElement("div");
            handle.id = "sidebar-resize-handle";
            handle.setAttribute("role", "separator");
            handle.setAttribute("aria-orientation", "vertical");
            handle.title = "Drag to resize sidebar";
            sideNav.appendChild(handle);
        }

        let startX = 0;
        let startWidth = 0;

        handle.addEventListener("pointerdown", function (event) {
            startX = event.clientX;
            startWidth = getCurrentWidth();

            handle.setPointerCapture(event.pointerId);
            document.body.classList.add("sidebar-resizing");

            event.preventDefault();
        });

        handle.addEventListener("pointermove", function (event) {
            if (!handle.hasPointerCapture(event.pointerId)) return;

            const delta = event.clientX - startX;
            applyWidth(startWidth + delta);
        });

        handle.addEventListener("pointerup", function (event) {
            if (handle.hasPointerCapture(event.pointerId)) {
              handle.releasePointerCapture(event.pointerId);
            }

            document.body.classList.remove("sidebar-resizing");
        });
    }

    const saved = parseInt(localStorage.getItem(storageKey), 10);
    if (!isNaN(saved)) {
       applyWidth(saved);
    }

    window.addEventListener("load", ensureResizeHandle);
})();
