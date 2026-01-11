# Makefile for VASP Syntax VS Code Extension

.PHONY: all install build package clean uninstall dev help

# Default target
all: build

# Install dependencies
node_modules: package.json
	npm install

# Build the extension (compile TypeScript)
build: node_modules
	npm run compile

# Build for development with watch mode
dev: node_modules
	npm run watch

# Package the extension as a .vsix file
package: build
	@command -v vsce >/dev/null 2>&1 || { echo "Installing vsce..."; npm install -g @vscode/vsce; }
	vsce package

# Install the extension to VS Code
install: package
	@echo "Installing VASP Syntax extension to VS Code..."
	@VSIX=$$(ls -t *.vsix 2>/dev/null | head -1); \
	if [ -z "$$VSIX" ]; then \
		echo "Error: No .vsix file found. Run 'make package' first."; \
		exit 1; \
	fi; \
	code --install-extension "$$VSIX" && echo "Extension installed successfully!"

# Install to VS Code Insiders
install-insiders: package
	@echo "Installing VASP Syntax extension to VS Code Insiders..."
	@VSIX=$$(ls -t *.vsix 2>/dev/null | head -1); \
	if [ -z "$$VSIX" ]; then \
		echo "Error: No .vsix file found. Run 'make package' first."; \
		exit 1; \
	fi; \
	code-insiders --install-extension "$$VSIX" && echo "Extension installed successfully!"

# Install by copying directly (alternative method without packaging)
install-dev: build
	@echo "Installing extension in development mode..."
	@EXTENSIONS_DIR="$$HOME/.vscode/extensions/vasp-syntax"; \
	mkdir -p "$$EXTENSIONS_DIR"; \
	cp -r package.json syntaxes language-configuration.json poscar-language-configuration.json out "$$EXTENSIONS_DIR/"; \
	echo "Extension installed to $$EXTENSIONS_DIR"; \
	echo "Restart VS Code to activate the extension."

# Uninstall the extension
uninstall:
	@echo "Uninstalling VASP Syntax extension..."
	code --uninstall-extension vasp-ase.vasp-syntax 2>/dev/null || true
	@EXTENSIONS_DIR="$$HOME/.vscode/extensions/vasp-syntax"; \
	if [ -d "$$EXTENSIONS_DIR" ]; then \
		rm -rf "$$EXTENSIONS_DIR"; \
		echo "Removed $$EXTENSIONS_DIR"; \
	fi
	@echo "Extension uninstalled. Restart VS Code."

# Clean build artifacts
clean:
	rm -rf out node_modules *.vsix

# Lint the TypeScript code
lint: node_modules
	npm run lint

# Run tests (if any)
test: build
	@echo "No tests configured yet"

# Show help
help:
	@echo "VASP Syntax VS Code Extension - Build Targets"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all              Build the extension (default)"
	@echo "  build            Compile TypeScript source"
	@echo "  dev              Build with watch mode for development"
	@echo "  package          Create .vsix package file"
	@echo "  install          Install extension to VS Code"
	@echo "  install-insiders Install extension to VS Code Insiders"
	@echo "  install-dev      Install directly without packaging (for development)"
	@echo "  uninstall        Remove extension from VS Code"
	@echo "  clean            Remove build artifacts"
	@echo "  lint             Run linter on TypeScript code"
	@echo "  test             Run tests"
	@echo "  help             Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - Node.js and npm"
	@echo "  - VS Code (for install targets)"
	@echo "  - vsce (installed automatically if needed)"
