// Extract each <svg> from the rationale page and render it to a high-DPI PNG
// with the light palette baked in, for embedding in the Word document.
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

const dir = __dirname;
const src = fs.readFileSync(path.join(dir, 'panel.html'), 'utf8');

// The token block, so the figures keep their colour meaning outside the page.
const tokens = src.match(/:root\{[\s\S]*?\}/)[0];
const svgs = src.match(/<svg\b[\s\S]*?<\/svg>/g);

(async () => {
  const browser = await chromium.launch({
    executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
  });
  const page = await browser.newPage({ deviceScaleFactor: 3 });

  for (let i = 0; i < svgs.length; i++) {
    const html = `<!doctype html><html><head><meta charset="utf-8">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
${tokens}
*{box-sizing:border-box}
body{margin:0;padding:14px;background:#ffffff;color:var(--ink);
     font-family:"IBM Plex Sans",system-ui,sans-serif;width:1180px}
svg{width:1152px;height:auto;display:block}
</style></head><body>${svgs[i]}</body></html>`;

    const file = path.join(dir, `fig-${i + 1}.html`);
    fs.writeFileSync(file, html);
    await page.goto('file://' + file);
    await page.waitForLoadState('networkidle').catch(() => {});
    await page.waitForTimeout(250);
    const el = await page.$('svg');
    await el.screenshot({ path: path.join(dir, `fig-${i + 1}.png`) });
    console.log(`fig-${i + 1}.png`);
  }

  await browser.close();
})();
