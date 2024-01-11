document.addEventListener('DOMContentLoaded', function() {
    // This will add the reset functionality to all the buttons with the id "submit".
    document.querySelectorAll('#submit').forEach(function(button) {
        button.addEventListener('click', function(event) {
            // Prevent the form from submitting to see the reset effect. 
            // Comment this out if you want the form to submit after resetting.
            event.preventDefault();

            // Reset all input fields.
            document.querySelectorAll('input[type="text"], input[type="email"]').forEach(function(input) {
                input.value = '';
            });

            // Reset all file inputs.
            document.querySelectorAll('input[type="file"]').forEach(function(input) {
                input.value = null;
            });

            // Reset all select dropdowns.
            document.querySelectorAll('select').forEach(function(select) {
                select.selectedIndex = 0;
            });

            // Empty the tissue list table.
            document.querySelector('#tBody').innerHTML = '';

            // Hide all alerts.
            document.querySelectorAll('.alert').forEach(function(alert) {
                alert.classList.add('d-none');
            });

            // Enable/Disable the necessary fields.
            document.querySelector('#global-file').disabled = true;
        });
    });
});
