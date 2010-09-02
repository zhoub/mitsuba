#include "previewsettingsdlg_cocoa.h"

@implementation CocoaRenderSettingsDlg
- (id) initWithParentView: (NSView *) view andDelegate: (PreviewSettingsDlg *) del {
	self = [super init];
	if (self) {
		parent = view;
		delegate = del;
		context = nil;
		if (![NSBundle loadNibNamed: @"PreviewSettings" owner: self])
			SLog(EError, "Unable to load PreviewSettings.nib");
		[panel setDelegate: self];
		active = NO;
	}
	return self;
}

- (void) hide {
	[panel orderOut: nil];
}

- (BOOL) isVisible {
	return [panel isVisible];
}

- (BOOL) isActiveWindow {
	return active;
}

- (void) showAt: (NSPoint) point {
	NSRect frame = [panel frame];
	CGFloat screenHeight = [[NSScreen mainScreen] visibleFrame].size.height;
	frame.origin.x = point.x - frame.size.width/2.0f;
	frame.origin.y = screenHeight - (point.y + frame.size.height/2.0f);
	[panel setFrame: frame display:YES];
	[panel makeKeyAndOrderFront: parent];
}

- (void) windowWillClose: (NSNotification *) notification {
	delegate->triggerClose();
}

- (void) windowDidBecomeKey: (NSNotification *) notification {
	active = YES;
}

- (void) windowDidResignKey: (NSNotification *) notification {
	active = NO;
}

- (void) updateUI {
	int sel_previewMethod = [previewMethod indexOfSelectedItem];
	int sel_toneMappingMethod = [toneMappingMethod indexOfSelectedItem];
	bool hasShadowMaps = (sel_previewMethod == EOpenGL || sel_previewMethod == EOpenGLSinglePass);
	bool hasKey = (sel_toneMappingMethod == EReinhard);
	bool hasGamma = [sRGB state] == NSOffState;
	bool hasDiffuseSources = [diffuseSourcesBox state] == NSOnState;
	std::string gammaStr = formatString("%.1f", [gamma floatValue]);
	std::string exposureStr = formatString("%.1f", [exposure floatValue]);

	[exposureLabel setStringValue: hasKey ? @"Burn :" : @"Exposure"];
	[shadowMapResolution setEnabled: (BOOL) hasShadowMaps];
	[[shadowMapLabel animator] setAlphaValue: hasShadowMaps ? 1.0f : 0.5f];
	[[shadowMapResolution animator] setAlphaValue: hasShadowMaps ? 1.0f : 0.5f];
	[reinhardKey setEnabled: (BOOL) hasKey];
	[[reinhardKeyLabel animator] setAlphaValue: hasKey ? 1.0f : 0.5f];
	[[reinhardKey animator] setAlphaValue: hasKey ? 1.0f : 0.5f];
	[gamma setEnabled: (BOOL) hasGamma];
	[[gammaLabel animator] setAlphaValue: hasGamma ? 1.0f : 0.5f];
	[[gamma animator] setAlphaValue: hasGamma ? 1.0f : 0.5f];
	[gammaValueLabel setStringValue: [NSString stringWithCString: gammaStr.c_str() encoding: NSASCIIStringEncoding]];
	[exposureValueLabel setStringValue: [NSString stringWithCString: exposureStr.c_str() encoding: NSASCIIStringEncoding]];
	[diffuseReceiversBox setEnabled: (BOOL) hasDiffuseSources];
	[[diffuseReceiversBox animator] setAlphaValue: hasDiffuseSources ? 1.0f : 0.5f];
}

- (IBAction) previewMethodChanged: (id) sender {
	delegate->triggerPreviewMethodChanged((EPreviewMethod) [previewMethod indexOfSelectedItem]);
	[self updateUI];
}

- (IBAction) toneMappingMethodChanged: (id) sender {
	EToneMappingMethod method = (EToneMappingMethod) [toneMappingMethod indexOfSelectedItem];
	if (method != context->toneMappingMethod) {
		delegate->triggerToneMappingMethodChanged(method);
		if (method == EGamma)
			[exposure setFloatValue: (float) context->exposure];
		else
			[exposure setFloatValue: (float) context->reinhardBurn];
	}
	[self updateUI];
}

- (IBAction) reinhardKeyChanged: (id) sender {
	delegate->triggerReinhardKeyChanged([reinhardKey floatValue]);
}

- (IBAction) exposureChanged: (id) sender {
	if ([toneMappingMethod indexOfSelectedItem] == EGamma)
		delegate->triggerExposureChanged([exposure floatValue]);
	else
		delegate->triggerReinhardBurnChanged([exposure floatValue]);
	[self updateUI];
}

- (IBAction) gammaChanged: (id) sender {
	delegate->triggerGammaChanged([sRGB state] == NSOnState, [gamma floatValue]);
	[self updateUI];
}

- (IBAction) sRGBChanged: (id) sender {
	delegate->triggerGammaChanged([sRGB state] == NSOnState, [gamma floatValue]);
	[self updateUI];
}

- (IBAction) diffuseSourcesChanged: (id) sender {
	delegate->triggerDiffuseSourcesChanged([diffuseSourcesBox state] == NSOnState);
	[self updateUI];
}

- (IBAction) diffuseReceiversChanged: (id) sender {
	delegate->triggerDiffuseReceiversChanged([diffuseReceiversBox state] == NSOnState);
}

- (IBAction) shadowMapResolutionChanged: (id) sender {
	int res, index = [shadowMapResolution indexOfSelectedItem];
	switch (index) {
		case 0: res = 64; break;
		case 1: res = 128; break;
		case 2: res = 256; break;
		case 3: res = 512; break;
		case 4: res = 1024; break;
		case 5: res = 2048; break;
		default:
			SLog(EError, "Invalid shadow map resolution!");
			return;
	}
	delegate->triggerShadowMapResolutionChanged(res);
}

- (IBAction) clampingChanged: (id) sender {
	delegate->triggerClampingChanged([clamping floatValue]);
}

- (IBAction) pathLengthChanged: (id) sender {
	delegate->triggerPathLengthChanged([pathLength intValue]);
}

- (IBAction) reset: (id) sender {
	[previewMethod selectItemAtIndex: EOpenGL];
	[toneMappingMethod selectItemAtIndex: EGamma];
	[shadowMapResolution selectItemAtIndex: 2];
	[pathLength setIntValue: 3];
	[clamping setFloatValue: 0.1];
	[gamma setFloatValue: 2.2];
	[exposure setFloatValue: 0.0];
	[reinhardKey setFloatValue: 0.18];
	[sRGB setState: NSOnState];
	[diffuseSourcesBox setState: NSOnState];
	[diffuseReceiversBox setState: NSOffState];
	[self previewMethodChanged: self];
	[self toneMappingMethodChanged: self];
	[self shadowMapResolutionChanged: self];
	[self pathLengthChanged: self];
	[self exposureChanged: self];
	[self clampingChanged: self];
	[self gammaChanged: self];
	[self diffuseSourcesChanged: self];
	[self diffuseReceiversChanged: self];
}

- (void) setContext: (SceneContext *) ctx {
	[previewMethod selectItemAtIndex: ctx->previewMethod];
	[toneMappingMethod selectItemAtIndex: ctx->toneMappingMethod];
	[reinhardKey setFloatValue: (float) ctx->reinhardKey];
	[pathLength setIntValue: ctx->pathLength];
	[clamping setFloatValue: (float) ctx->clamping];
	[gamma setFloatValue: (float) ctx->gamma];
	[sRGB setState: ctx->srgb ? NSOnState : NSOffState];
	[diffuseSourcesBox setState: ctx->diffuseSources ? NSOnState : NSOffState];
	[diffuseReceiversBox setState: ctx->diffuseReceivers ? NSOnState : NSOffState];

	if (ctx->toneMappingMethod == EGamma)
		[exposure setFloatValue: (float) ctx->exposure];
	else
		[exposure setFloatValue: (float) ctx->reinhardBurn];

	int shadowMapIdx;
	switch (ctx->shadowMapResolution) {
		case 64: shadowMapIdx = 0; break;
		case 128: shadowMapIdx = 1; break;
		case 256: shadowMapIdx = 2; break;
		case 512: shadowMapIdx = 3; break;
		case 1024: shadowMapIdx = 4; break;
		case 2048: shadowMapIdx = 5; break;
		default:
			SLog(EError, "Invalid shadow map resolution!");
			return;
	}

	[shadowMapResolution selectItemAtIndex: shadowMapIdx];
	[self updateUI];
	context = ctx;
}

@end

PreviewSettingsDlg::PreviewSettingsDlg(QWidget *parent) : QObject(parent) {
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	m_dlg = [[CocoaRenderSettingsDlg alloc] 
		initWithParentView: (NSView *) parent->winId()
		andDelegate: this
	];
	[pool release];
}

PreviewSettingsDlg::~PreviewSettingsDlg() {
}

void PreviewSettingsDlg::show() {
	QWidget *widget = (QWidget *) parent();
	QPoint pos = widget->pos();
	QSize size = widget->size();
	[m_dlg showAt: NSMakePoint(pos.x() + size.width()/2, 
		pos.y() + size.height()/2)];
}

void PreviewSettingsDlg::hide() {
	[m_dlg hide];
}

bool PreviewSettingsDlg::isVisible() {
	return (bool) [m_dlg isVisible];
}

void PreviewSettingsDlg::setContext(SceneContext *ctx) {
	[m_dlg setContext: ctx];
}

bool PreviewSettingsDlg::isActiveWindow() const {
	return [m_dlg isActiveWindow];
}
