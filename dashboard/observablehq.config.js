export default {
  title: "C. elegans TE Explorer",
  root: "src",
  base: "/Low-input-ribosome-profiling-of-Celegans-embryogenesis/",
  theme: "air",
  sidebar: true,
  toc: false,
  pages: [
    {name: "Explorer", path: "/"},
    {name: "About", path: "/about"}
  ],
  head: `<style>
    :root {
      --theme-background: #fafafa;
      --accent: #D64545;
      --ribo: #E0662B;
      --rna: #3A7CA5;
      --ref-grey: #9AA0A6;
      --ink: #1A2027;
      --rule: #E6E8EB;
    }
    #observablehq-main {
      color: var(--ink);
    }
    #observablehq-main h1,
    #observablehq-main h2,
    #observablehq-main h3 {
      color: var(--ink);
      letter-spacing: -0.01em;
      line-height: 1.15;
    }
    .card {
      background: white;
      border: 1px solid #eaeaea;
      box-shadow: 0 2px 6px rgba(0,0,0,0.04);
      border-radius: 8px;
    }
    .card p {
      max-width: 100%;
    }
    .eyebrow {
      font-size: 11px;
      font-weight: 600;
      letter-spacing: 0.09em;
      text-transform: uppercase;
      color: #8A9098;
      margin: 0 0 4px;
    }
    #observablehq-main p,
    #observablehq-main h1,
    #observablehq-main h2,
    #observablehq-main h3,
    #observablehq-main table,
    #observablehq-main blockquote,
    #observablehq-main ul,
    #observablehq-main ol,
    #observablehq-main hr {
      max-width: 100%;
    }
  </style>`
};
