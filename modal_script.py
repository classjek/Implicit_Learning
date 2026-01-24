import modal

app = modal.App("test-culorads")

# Create image with CuLoRADS binary
cuda_image = (
    modal.Image.from_registry("nvidia/cuda:12.6.1-base-ubuntu22.04")
    .apt_install("wget", "tar")
    .run_commands(
        # Download binary
        "wget -P /opt https://github.com/COPT-Public/cuLoRADS/releases/download/v1.0.0/cuLoRADS.tar.gz",
        # Extract it
        "cd /opt && tar -xzvf cuLoRADS.tar.gz",
        # CHeck executability
        "chmod +x /opt/cuLoRADS/bin/cuLoRADS",
    )
)

@app.function(gpu="A10G", image=cuda_image)
def check_culorads():
    """Check if CuLoRADS is installed and working"""
    import subprocess
    
    # check binary
    result = subprocess.run(
        ["ls", "-lh", "/opt/cuLoRADS/bin/cuLoRADS"],
        capture_output=True,
        text=True
    )
    print("Binary location:")
    print(result.stdout)
    
    # confirm running the binary
    result = subprocess.run(
        ["/opt/cuLoRADS/bin/cuLoRADS"],
        capture_output=True,
        text=True,
        timeout=30
    )
    print("\nCuLoRADS output:")
    print(result.stdout)
    print(result.stderr)
    
    return "Check complete"

if __name__ == "__main__":
    with app.run():
        check_culorads.remote()