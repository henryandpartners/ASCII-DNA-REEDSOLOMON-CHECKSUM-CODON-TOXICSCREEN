document.addEventListener('DOMContentLoaded', () => {
    const inputText = document.getElementById('input-text');
    const outputText = document.getElementById('output-text');
    const encodeBtn = document.getElementById('encode-btn');
    const decodeBtn = document.getElementById('decode-btn');
    const paddingInput = document.getElementById('padding-input');
    const imageUpload = document.getElementById('image-upload');
    const convertBtn = document.getElementById('convert-btn');
    const charCount = document.getElementById('char-count');
    const basePairCount = document.getElementById('base-pair-count');

    // --- Event Listeners ---

    inputText.addEventListener('input', () => {
        charCount.textContent = inputText.value.length;
    });

    encodeBtn.addEventListener('click', () => {
        const text = inputText.value;
        fetch('/encode', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ text: text })
        })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                outputText.value = data.dna;
                paddingInput.value = data.padding;
                basePairCount.textContent = data.base_pairs;
            } else {
                outputText.value = `Error: ${data.error}`;
            }
        });
    });

    decodeBtn.addEventListener('click', () => {
        const dna = inputText.value;
        const padding = paddingInput.value;
        fetch('/decode', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ dna: dna, padding: padding })
        })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                outputText.value = data.text;
                basePairCount.textContent = 0;
            } else {
                outputText.value = `Error: ${data.error}`;
            }
        });
    });

    convertBtn.addEventListener('click', () => {
        const file = imageUpload.files[0];
        if (!file) {
            inputText.value = 'Please select an image file first.';
            return;
        }

        const formData = new FormData();
        formData.append('image', file);

        fetch('/image-to-ascii', {
            method: 'POST',
            body: formData
        })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                inputText.value = data.ascii;
                charCount.textContent = inputText.value.length;
            } else {
                inputText.value = `Error: ${data.error}`;
            }
        });
    });
});
