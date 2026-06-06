// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
/**
Set the default toggle level to 1
*/

function toggleLevel(level) {
    document.querySelectorAll('tr[id^="row_"]').forEach(function(row) {
        const depth = row.id.replace(/^row_/, '').split('_').filter(Boolean).length - 1;
        row.style.display = depth < level ? '' : 'none';
    });
}

document.addEventListener('DOMContentLoaded', function() {
    toggleLevel(1);
});
