import { mutation, query } from "./_generated/server";
import { v } from "convex/values";

const structureRecordArgs = {
  tenant_id: v.string(),
  workspace_id: v.string(),
  structure_id: v.string(),
  source: v.string(),
  cif: v.string(),
  structure: v.any(),
  params: v.optional(v.any()),
  created_at: v.string(),
  updated_at: v.string(),
  raw_input_chunk_count: v.number(),
} as const;

export const putStructureRecord = mutation({
  args: structureRecordArgs,
  handler: async (ctx, args) => {
    const existing = await ctx.db
      .query("structures")
      .withIndex("by_key", (q) =>
        q
          .eq("tenant_id", args.tenant_id)
          .eq("workspace_id", args.workspace_id)
          .eq("structure_id", args.structure_id),
      )
      .first();

    if (existing) {
      await ctx.db.patch(existing._id, {
        source: args.source,
        cif: args.cif,
        structure: args.structure,
        params: args.params,
        updated_at: args.updated_at,
        raw_input_chunk_count: args.raw_input_chunk_count,
      });
      return { updated: true } as const;
    }

    await ctx.db.insert("structures", args);
    return { updated: false } as const;
  },
});

export const putStructureRawChunk = mutation({
  args: {
    tenant_id: v.string(),
    workspace_id: v.string(),
    structure_id: v.string(),
    chunk_index: v.number(),
    chunk_text: v.string(),
  },
  handler: async (ctx, args) => {
    const existing = await ctx.db
      .query("structure_raw_chunks")
      .withIndex("by_structure", (q) =>
        q
          .eq("tenant_id", args.tenant_id)
          .eq("workspace_id", args.workspace_id)
          .eq("structure_id", args.structure_id)
          .eq("chunk_index", args.chunk_index),
      )
      .first();

    if (existing) {
      await ctx.db.patch(existing._id, { chunk_text: args.chunk_text });
      return { updated: true } as const;
    }

    await ctx.db.insert("structure_raw_chunks", args);
    return { updated: false } as const;
  },
});

export const getStructureRecord = query({
  args: {
    tenant_id: v.string(),
    workspace_id: v.string(),
    structure_id: v.string(),
  },
  handler: async (ctx, args) => {
    return ctx.db
      .query("structures")
      .withIndex("by_key", (q) =>
        q
          .eq("tenant_id", args.tenant_id)
          .eq("workspace_id", args.workspace_id)
          .eq("structure_id", args.structure_id),
      )
      .first();
  },
});

export const getStructureRecordById = query({
  args: {
    structure_id: v.string(),
  },
  handler: async (ctx, args) => {
    return ctx.db
      .query("structures")
      .withIndex("by_structure_id", (q) =>
        q.eq("structure_id", args.structure_id),
      )
      .first();
  },
});

export const getStructureRawChunks = query({
  args: {
    tenant_id: v.string(),
    workspace_id: v.string(),
    structure_id: v.string(),
  },
  handler: async (ctx, args) => {
    return ctx.db
      .query("structure_raw_chunks")
      .withIndex("by_structure", (q) =>
        q
          .eq("tenant_id", args.tenant_id)
          .eq("workspace_id", args.workspace_id)
          .eq("structure_id", args.structure_id),
      )
      .collect();
  },
});
